# -*- coding: utf-8 -*-
"""A module to read xyz-formatted files.

This module reads files in .xyz format.  This file can also be run as a script to convert xyz files
into several other formats.  Currently supported formats are csv and cel.  Written by Charles
Evans.

Todo:
    * Implement other output formats (cif, maybe z-matrix)
    * Command line argument version w/o live user interaction
"""


class XYZ:
    """A class to contain the data of an .xyz file in an accessible format.

    Attributes:
        num_sites (int): The number of described atomic sites
        comment_string (str): An arbitraty comment string
        sites (list((E,x,y,z))): A list of tuples containing site information (E is the element)
    """

    def __init__(self, num_sites, comment_string, sites):
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
    """Exports an xyz-formatted file.

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


####################
# Script Mode Only #
####################

def define_file(in_or_out, file_ext, inherit_name=""):
    """Handles interaction with user to get a valid filename for reading or writing.

    Args:
        in_or_out (str): Must be one of "input" or "output"
        file_ext (str): File extension string, including leading period (e.g. ".csv" not "csv")
        inherit_name (str, optional): Default file name if user provides none

    Returns:
        filename (str): An OS-valid path to a file with the extension defined by file_ext

    Todo:
        * Implement checks for args, throw errors if they're malformatted
    """
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


def progressBar(iterable, prefix='', suffix='', decimals=1, length=50, fill='█', printEnd="\r"):
    """Call in a loop to create terminal progress bar.

       Written by stackoverflow user Greenstick.

    Args:
        iterable (Iterable): Iterable object
        prefix (str, optional): Prefix string
        suffix (str, optional): Suffix string
        decimals (int, optional): Positive number of decimals in percent complete
        length (int, optional): Character length of bar
        fill (str, optional): Bar fill character
        printEnd (str, optional): End character (e.g. "\r", "\r\n")
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


def main():
    """Interactive script mode.

    Todo:
        * Move the different file type operations into their own methods
            (for readability & maintainability)
    """
    while True:
        infile = define_file("input", ".xyz")
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

    # This list of supported output formats needs to be kept up to date
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
        if out_type == "exit":
            sys.exit()
        if out_type not in supported_formats:
            print("Not a recognized format")
            continue
        # All other cases should be supported formats

        #######
        # CSV #
        #######
        # See RFC 4180 for a description of the csv format
        if out_type == "csv":
            from csv import writer
            print("Writing xyz data to csv file...")
            outfile = define_file("output", ".csv", inherit_name=infile)
            delim = str(input("Please choose a delimiting character (comma default): "))
            if delim == "":
                delim = ","  # Pressing enter with no other input returns the default value
            with open(outfile, 'w', newline='') as file:
                csvwriter = writer(file, delimiter=delim)  # This could fail for weird delimiters
                csvwriter.writerow(["Element", "x", "y", "z"])
                for site in progressBar(xyz_obj.sites, prefix="Writing:"):
                    csvwriter.writerow(site)

        #######
        # CEL #
        #######
        # See https://er-c.org/barthel/drprobe/celfile.html for a description of the cel format
        if out_type == "cel":
            print("Writing data to cel file...")
            outfile = define_file("output", ".cel", inherit_name=infile)
            if str(input("Would you like to define Debey-Waller paramters?"
                         " (y/n)\n")).lower() in ["y", "yes"]:
                deb_switch = True
            else:
                deb_switch = False

            cel_comment = "Number of atoms: " + str(xyz_obj.num_sites) \
                + "; " + xyz_obj.comment_string + "\n"
            deb_dict = {}  # We only need this if deb_wal == True, but it's easier to just make it
            a_max, a_min, b_max, b_min, c_max, c_min = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            alpha, beta, gamma = 90, 90, 90  # Not sure if non-ortho coordinates are possible
            for site in progressBar(xyz_obj.sites, prefix="Calculating span:"):
                if site[0] not in deb_dict:
                    deb_dict[site[0]] = 0.0  # Default D-W param is 0
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

            if deb_switch:
                for elem, _ in deb_dict.items():
                    while True:
                        try:
                            dw = float(input("Please input the D-W parameter for " +
                                             elem + " in A^2 (or 0 to skip): ")) / 100
                            if dw < 0:
                                print("Please enter a value greater than or equal to zero.")
                                continue
                            deb_dict.update({elem: round(dw, 6)})
                            break
                        except ValueError:
                            print("Please enter a numeric value.")
                            continue

            with open(outfile, 'w') as file:
                file.write(cel_comment)  # The first line in a cel file is a comment string
                header = "0 " + str(a_span) + " " + str(b_span) + " " + str(c_span) + " " \
                    + str(alpha) + " " + str(beta) + " " + str(gamma) + "\n"
                file.write(header)  # The second line in a cel file defines its extent

                for site in progressBar(xyz_obj.sites, prefix="Writing"):
                    x = round((float(site[1])/10 + abs(a_min)) / a_span, 6)
                    y = round((float(site[2])/10 + abs(b_min)) / b_span, 6)
                    z = round((float(site[3])/10 + abs(c_min)) / c_span, 6)
                    line = site[0] + " " + str(x) + " " + str(y) + " " + str(z) \
                        + " 1.0 " + str(deb_dict[site[0]]) + " 0.0 0.0 0.0\n"
                    file.write(line)


if __name__ == "__main__":
    print("XYZReader running as interactive script...")
    import os
    import sys
    main()
