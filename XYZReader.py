# -*- coding: utf-8 -*-
"""A module to read xyz-formatted files.

This module reads files in .xyz format.  This file can also be run as a script to convert xyz files
into several other formats.  Currently supported formats are csv and cel.  Written by Charles
Evans.

Todo:
    * Import from other formats and output to xyz
    * Implement other output formats (cif, maybe z-matrix)
    * Command line argument version w/o live user interaction
"""


class XYZ:
    """A class to contain the data of an .xyz file in an accessible format.

    Attributes:
        num_sites (int): The number of described atomic sites

        comment_string (str): An arbitraty comment string

        sites (list(E,x,y,z)): A list of tuples containing site information (E is the element)
    """

    def __init__(self, num_sites, comment_string, sites):
        """Initialize an XYZ object."""
        self.num_sites = num_sites
        self.comment_string = comment_string
        self.sites = sites


########################
# Importable Functions #
########################


def readXYZ(file):
    """Read in an .xyz formatted file, return an instance of XYZ class.

    Args:
        file (str): The path to the xyz file to be read

    Returns:
        XYZ: An initialized XYZ object containing the information from the file
    """
    lines = file.readlines()
    # .xyz files have two lines before the atomic positions are described:
    #   - the first line is the number of atoms described in the file
    #   - the second line is a comment string
    num_sites = int(lines.pop(0))
    comment_string = lines.pop(0)
    # The test of the lines in the document contain identities and locations
    #  of atoms
    sites = []
    for line in lines:
        cols = line.split()
        # By default the x, y, z coords will be strings, but we want floats
        for i in range(1, 4):
            cols[i] = float(cols[i])
        site = tuple(cols)
        sites.append(site)
    return XYZ(num_sites, comment_string, sites)


def writeXYZ(xyz_obj, outfile):
    """Export an xyz-formatted file.

    Args:
        xyz_obj (XYZ): An initialized XYZ object containing the information to be written

        outfile (str): The path to the file which will be written
    """
    with open(outfile, 'w') as file:
        file.write(xyz_obj.num_sites, "\n")
        file.write(xyz_obj.comment_string, "\n")
        for site in xyz_obj.sites:  # No prog bar in importable, do it in the script if desired
            line = " ".join(site) + "\n"
            file.write(line)


##################
# Live Mode Only #
##################

def _interactive_define_file(in_or_out, file_ext, inherit_name=""):
    # Interactive-mode file definition
    # TODO: Implement checks for args, throw errors if they're malformatted
    while True:
        if in_or_out == "input":
            filename = str(input("Please enter the input file location: "))
            if not filename.endswith(file_ext):
                filename = filename + file_ext
            if not os.path.isfile(filename):
                if str(input("File does not exist. Try again? (y/n)\n")).lower() in ["y", "yes"]:
                    continue
                else:
                    sys.exit()
            else:
                return filename

        if in_or_out == "output":
            while True:
                filename = input("Please name your output file (blank for same name): ").strip()
                if filename == "":  # Use the same filename as the currently loaded xyz file
                    # For this to work, must pass optional param inherit_name with input filename
                    filename = inherit_name[:-4]  # Trim extension from inherited name
                if not filename.endswith(file_ext):
                    filename = filename + file_ext
                if os.path.isfile(filename):
                    if str(input("Existing file; overwrite? (y/n)\n")).lower() in ["y", "yes"]:
                        return filename
                    elif str(input("Try again? (y/n)\n")).lower() in ["y", "yes"]:
                        continue
                    else:
                        sys.exit()
                return filename


def _parse_args(args):
    parser = argparse.ArgumentParser(
        description="Convert back and forth from the .xys file format")
    parser.add_argument("--interactive", action="store_true",
                        help="launch in interactive mode, ignoring other options")
    parser.add_argument("infile", nargs="?", default="",
                        help="the path to the input file; currently only supports xyz files")
    parser.add_argument("outfile", nargs="?", default="",
                        help="the path to the output file")
    parser.add_argument("--type", action="append", nargs="+", choices=["cel", "csv"],
                        help="the type of file to output; multiple types may be listed at once; defaults to inferring from file extension of outfile")
    parser.add_argument("--overwrite", action="store_true",
                        help="allow outfile to be overwritten if it already exists")
    # Add additional filetypes to choices as they're implemented
    opts = parser.parse_args()
    return opts


def _progressBar(iterable, prefix='', suffix='', decimals=1, length=50, fill='█', printEnd="\r"):
    r"""Call in a loop to create terminal progress bar.

       Written by stackoverflow user Greenstick.

    Args:
        iterable (Iterable): Iterable object

        prefix (str, optional): Prefix string

        suffix (str, optional): Suffix string

        decimals (int, optional): Positive number of decimals in percent complete

        length (int, optional): Character length of bar

        fill (str, optional): Bar fill character

        printEnd (str, optional): End character (e.g. "\r", "\r \n")
    """
    total = len(iterable)

    # Progress Bar Printing Function
    def printProgressBar(iteration):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=printEnd)
    # Initial Call
    print()  # Newline before progress bar
    printProgressBar(0)
    # Update Progress Bar
    for i, item in enumerate(iterable):
        yield item
        printProgressBar(i + 1)
    # Print newline and completion message
    print()
    print("Complete!")


def _write_csv(outfile, xyz_obj, delim=","):
    # See RFC 4180 for a description of the csv format
    from csv import writer
    print("Writing xyz data to csv file...")
    with open(outfile, 'w', newline='') as file:
        csvwriter = writer(file, delimiter=delim)  # This could fail for weird delimiters
        csvwriter.writerow(["Element", "x", "y", "z"])
        for site in _progressBar(xyz_obj.sites, prefix="Writing:"):
            csvwriter.writerow(site)


def _write_cel(outfile, xyz_obj, deb_dict={}):
    # See https://er-c.org/barthel/drprobe/celfile.html for a description of the cel format
    print("Writing data to cel file...")
    cel_comment = "Number of atoms: " + str(xyz_obj.num_sites) \
        + "; " + xyz_obj.comment_string + "\n"
    a_max, a_min, b_max, b_min, c_max, c_min = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    alpha, beta, gamma = 90, 90, 90  # Not sure if non-ortho coordinates are possible
    # TODO: test this and find out?
    for site in _progressBar(xyz_obj.sites, prefix="Calculating span:"):
        site.append(deb_dict.get(site[0], 0.0))  # Default D-W param is 0
        # TODO: there must be a more elegant way to write this part...
        if site[1] > a_max:
            a_max = site[1]
        if site[1] < a_min:
            a_min = site[1]
        if site[2] > b_max:
            b_max = site[2]
        if site[2] < b_min:
            b_min = site[2]
        if site[3] > c_max:
            c_max = site[3]
        if site[3] < c_min:
            c_min = site[3]
    a_span, b_span, c_span = a_max-a_min, b_max-b_min, c_max-c_min
    # xyz files are in A, but we need spans in nm for cel files
    a_span, b_span, c_span = round(a_span/10, 6), round(b_span/10, 6), round(c_span/10, 6)

    # TODO: copy this block into the interactive section of _main
    # if deb_switch:
    #     for elem, _ in deb_dict.items():
    #         while True:
    #             try:
    #                 dw = float(input("Please input the D-W parameter for " +
    #                                  elem + " in A^2 (or 0 to skip): ")) / 100
    #                 if dw < 0:
    #                     print("Please enter a value greater than or equal to zero.")
    #                     continue
    #                 deb_dict.update({elem: round(dw, 6)})
    #                 break
    #             except ValueError:
    #                 print("Please enter a numeric value.")
    #                 continue

    with open(outfile, 'w') as file:
        file.write(cel_comment)  # The first line in a cel file is a comment string
        header = "0 " + str(a_span) + " " + str(b_span) + " " + str(c_span) + " " \
            + str(alpha) + " " + str(beta) + " " + str(gamma) + "\n"
        file.write(header)  # The second line in a cel file defines its extent

        for site in _progressBar(xyz_obj.sites, prefix="Writing"):
            x = round((float(site[1])/10 + abs(a_min)) / a_span, 6)
            y = round((float(site[2])/10 + abs(b_min)) / b_span, 6)
            z = round((float(site[3])/10 + abs(c_min)) / c_span, 6)
            line = site[0] + " " + str(x) + " " + str(y) + " " + str(z) \
                + " 1.0 " + str(site[-1]) + " 0.0 0.0 0.0\n"
            file.write(line)


def _main(opts):
    # TODO: Move the different file type operations into their own methods
    # TODO: Full rewrite to support opts!

    ####################
    # Interactive Mode #
    ####################
    if opts.interactive:
        while True:
            infile = _interactive_define_file("input", ".xyz")
            try:
                with open(infile) as file:
                    xyz_obj = readXYZ(file)
            except Exception as ex:
                print("Something went wrong while reading from ", infile)
                sys.exit(ex.message)
            break
        # Print some info about the file so the user can confirm it's the right one
        print("File ", infile, ": ", xyz_obj.comment_string)
        print("File contains ", xyz_obj.num_sites, " atomic sites")
        _ = input("Press Enter to confirm...")
        supported_formats = {
            "exit": "Exit the script",
            "cel": "Dr Probe super-cell file",
            "csv": "Generic comma-separated value file"}
        while True:
            out_type = str(input("Please select the file type to export to (h for help): "))
            # Handle special cases: help request, not a supported format
            if out_type in ["h", "help"]:
                print("Supported filetypes:\n")
                for key, val in supported_formats.items():
                    print(key, " : ", val)
                continue
            elif out_type == "exit":
                sys.exit()
            elif out_type not in supported_formats:
                print("Not a recognized format")
                continue
            # All other cases should be supported formats
            if "csv" in out_type:
                pass

            if "cel" in out_type:
                pass

    #################
    # Argument Mode #
    #################
    else:
        infile = opts.inflie
        # Check that the infile exists
        if not infile.endswith(".xyz"):
            infile = infile + ".xyz"
        if not os.path.isfile(infile):
            print("File " + infile + " does not exist")
            sys.exit()
        # Create the xyz_obj from the data in infile
        try:
            with open(infile) as file:
                xyz_obj = readXYZ(file)
        except Exception as ex:
            print("Something went wrong while reading from ", infile)
            sys.exit(ex.message)

        # Safety check outfile (no overwriting without permission!)
        outfile = opts.outfile
        if os.path.isfile(outfile) and not opts.overwrite:
            print("File " + outfile +
                  " already exists. Pass the --overwrite option to overwrite it.")
            sys.exit()
        # Split off and ext; ext handled by out_type
        outfile_name, outfile_ext = os.path.splitext(outfile)
        out_type = opts.type
        # Default to inferring form outfile extension
        if out_type == []:
            if "csv" in outfile_ext:
                out_type = ["csv"]
            if "cel" in outfile_ext:
                out_type = ["cel"]

        if "csv" in out_type:
            _write_csv(outfile_name+".csv", xyz_obj)

        if "cel" in out_type:
            _write_cel(outfile_name+".cel", xyz_obj)


if __name__ == "__main__":
    import os
    import sys
    import argparse
    opts = _parse_args(sys.argv)
    # _main(opts)
