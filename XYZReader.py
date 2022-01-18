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

####################
# Importable Class #
####################

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
    # TODO: Sanity checking (file format is as expected) should be implemented
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


#################
# CMD Mode Only #
#################

def _parse_args(args):
    parser = argparse.ArgumentParser(
        description="Convert back and forth from the .xyz file format")
    parser.add_argument("infile", nargs="+", default="",
                        help="The path to the input file; currently only supports xyz files.  Multiple files can be passed at once for batch conversion.")
    parser.add_argument("--outfile", nargs="?", default="",
                        help="The path to the output file.  If not specified, uses the same path as the input file with a new extension appropriate to --type.")
    parser.add_argument("--type", action="append", nargs="+", choices=["cel", "csv"],
                        help="The type of file to output; multiple types may be listed at once, or type can be inferred from file extension of --outfile if provided.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Allow file overwriting on output if the output file already exists (default is to not allow overwriting).")
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


def _write_csv(outfile, xyz_obj, overwrite, delim=","):
    # See RFC 4180 for a description of the csv format
    from csv import writer
    
    # Safety check outfile (no overwriting without permission!)
    if os.path.isfile(outfile) and not overwrite:
        print("File " + outfile +
              " already exists. Use the --overwrite option if you want to overwrite it.")
        sys.exit()

    print("Writing xyz data to csv file...")
    with open(outfile, 'w', newline='') as file:
        csvwriter = writer(file, delimiter=delim)  # This could fail for weird delimiters
        csvwriter.writerow(["Element", "x", "y", "z"])
        for site in _progressBar(xyz_obj.sites, prefix="Writing:"):
            csvwriter.writerow(site)


def _write_cel(outfile, xyz_obj, overwrite, deb_dict={}):
    # See https://er-c.org/barthel/drprobe/celfile.html for a description of the cel format
    
    # Safety check outfile (no overwriting without permission!)
    if os.path.isfile(outfile) and not overwrite:
        print("File " + outfile +
              " already exists. Use the --overwrite option if you want to overwrite it.")
        sys.exit()

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
    infile = opts.infile[0]
    # TODO: When non-xyz inputs are supported, this part needs to change
    # Check that the infile exists (append .xyz if needed) and read in data
    if not infile.endswith(".xyz"):
        infile = infile + ".xyz"
    if not os.path.isfile(infile):
        print("File " + infile + " does not exist.")
        sys.exit()
    # Create the xyz_obj from the data in infile
    try:
        with open(infile) as file:
            xyz_obj = readXYZ(file)
    except Exception as ex:
        print("Something went wrong while reading from ", infile)
        sys.exit(ex.message)

    # CASE: Outfile was specified as argument
    if not opts.outfile == []:
        outfile = opts.outfile
        # Split off ext; ext handled by out_type
        outfile_name, outfile_ext = os.path.splitext(outfile)
        if not opts.type == []:
            out_type = opts.type
        # Default to inferring form outfile extension if type not passed
        else:
            out_type = outfile_ext.replace(".", "")

    # CASE: Outfile is not specified; use the same path & name as infile but with new ext
    else:
        outfile, _ = os.path.splitext(infile)
        if not opts.type == []:
            out_type = opts.type
        # In this scenario, type must be passed.  If it isn't, error and exit.
        else:
            print("Output type could not be inferred.  Please specify --type or explicitly define --outfile.")
            sys.exit()
    
    # Set out_type
    #TODO: There's probably a more elegant way to implement this that will be easier to upkeep
    if out_type == []:
        if "csv" in outfile_ext:
            out_type = ["csv"]
        if "cel" in outfile_ext:
            out_type = ["cel"]

    # Write out the data in the correct format
    # TODO: There's probably a more elegant way to do this too; combine with above?
    if "csv" in out_type:
        _write_csv(outfile_name+".csv", xyz_obj, opts.overwrite)

    if "cel" in out_type:
        _write_cel(outfile_name+".cel", xyz_obj, opts.overwrite)


# Run main loop iff called from the command line, but not when imported as module
if __name__ == "__main__":
    import os
    import sys
    import argparse
    opts = _parse_args(sys.argv)
    _main(opts)
