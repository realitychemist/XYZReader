"""
This module reads files in .xyz format.  This file can also be run as a script to convert xyz files
into several other formats.  Currently supported export formats are xyz, csv, and cel

:Todo:
    * Test this! I don't know if everything is working right now.
    * When CrystalMaker exports XYZ from a structure made of multiple models, each gets its own
        number of atoms and comment string (blank); must recognize and handle this situation
    * Implement other output formats (cif, FEFF ATOMS card)
"""
import typing
import os
if typing.TYPE_CHECKING:
    from argparse import Namespace

####################
# Importable Class #
####################


class XYZ:
    """
    A class to contain the data of an .xyz file in an accessible format.

    :Attributes:
        :num_sites: Number of atomic sites represented in the object (should be equal to len(sites))
        :comment_string: An arbitraty comment string
        :sites: Contians site information tuples (E, x, y, z) where E is the element shortname
            and x, y, and z are coordinate positions in angstroms
    """

    def __init__(self, num_sites: int, comment_string: str, sites: list):
        """Initialize an XYZ object."""
        self.num_sites = num_sites
        self.comment_string = comment_string
        self.sites = sites


########################
# Importable Functions #
########################


def read_xyz(file: typing.TextIO, symbol_first: bool = True) -> XYZ:
    """
    Reads in an .xyz formatted file, return an instance of XYZ class.
    :Params:
        :file: An xyz-formatted file object, for example the result of open(/path/to/file.xyz)
        :symbol_first (bool): Optional, sets whether the element symbols are listed first (True) or last (False)
            for each site (default is True)
    :Returns:
        An initialized XYZ object containing the information from the file
    """
    lines = file.readlines()
    # .xyz files have two lines before the atomic positions are described:
    #   - the first line is the number of atoms described in the file
    #   - the second line is a comment string
    # Now with sanity checking!

    err_msg = {"l1": "Malformatted XYZ: first line should contain number of sites, which must be a positive integer",
               "sval": "Malformatted XYZ: file contains site coordinates which are not castable to floating point",
               "sid": "Site contains an unrecognized element symbol"}
    known_elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
                      "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
                      "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
                      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
                      "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
                      "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
                      "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
                      "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
                      "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
    try:
        num_sites = int(lines.pop(0))
        if num_sites <= 0:
            raise ValueError(err_msg["l1"])
        if len(lines) < num_sites:
            raise UserWarning("Fewer sites than expected, something might be wrong with the file")
        if len(lines) > num_sites:
            raise UserWarning("More sites than expected, something might be wrong with the file")
    except ValueError:
        print(err_msg["l1"])

    comment_string = lines.pop(0)
    # Probably the most useful thing we can do with the comment string is print it, to help users
    # tell if something has gone wrong
    print(f"Reading XYZ file with comment string: {comment_string}")

    # The rest of the lines in the document contain identities and locations of atoms
    sites = []
    for ln, line in enumerate(lines):
        cols = line.split()

        if symbol_first:
            if cols[0] not in known_elements:
                raise UserWarning(err_msg["sid"] + f" on line {ln}: " + line)
        else:  # symbols last
            if cols[-1] not in known_elements:
                raise UserWarning(err_msg["sid"] + f" on line {ln}: " + line)

        try:
            if symbol_first:
                for i in range(1, 4):
                    cols[i] = float(cols[i])
            else:  # symbols last
                for i in range(0, 3):
                    cols[i] = float(cols[i])
        except TypeError:
            print(err_msg["sval"] + f" on line {ln}: " + line)

        site = tuple(cols)
        sites.append(site)
    # noinspection PyUnboundLocalVariable
    return XYZ(num_sites, comment_string, sites)


def write_xyz(xyz_obj: XYZ, outfile: str | os.PathLike, overwrite: bool = False):
    """
    Exports an xyz-formatted file.

    :Params:
        :xyz_obj: An initialized XYZ object containing the information to be written
        :outfile: The path to the file which will be written, including the file extension
    """
    _overwrite_check(outfile, overwrite)
    with open(outfile, 'w') as file:
        file.write(f"{xyz_obj.num_sites}\n")
        file.write(f"{xyz_obj.comment_string}\n")
        for site in xyz_obj.sites:
            line = " ".join(site) + "\n"
            file.write(line)


def write_csv(xyz_obj: XYZ, outfile: str | os.PathLike, overwrite: bool = False, **kwargs):
    """
    Writes out a CSV file based on the data in an XYZ object. See RFC 4180 for a description of
    the csv format.

    :Params:
        :xyz_obj: Contains the data to be written
        :outfile: The path to the outfile, including extension
        :overwrite: Optional, disables overwriting if False (default is False)
        :**kwargs: Passed to csv.writer
    """
    _overwrite_check(outfile, overwrite)
    from csv import writer

    with open(outfile, 'w', newline='') as file:
        csvwriter = writer(file, **kwargs)
        csvwriter.writerow(["Element", "x", "y", "z"])
        for site in xyz_obj.sites:
            csvwriter.writerow(site)


def write_cel(xyz_obj: XYZ, outfile: str | os.PathLike,
              uc_params: tuple[float, float, float], multiplicity: tuple[int, int, int] = (1, 1, 1),
              overwrite: bool = False, deb_dict: dict[str, float] | None = None):
    """
    Writes out a CEL file based on the data in an XYZ object.  See https://er-c.org/barthel/drprobe/celfile.html
    for a description of the cel format.

    :Params:
        :xyz_obj: Contains the data to be written
        :outfile: The path to the outfile, including extension
        :uc_params: Orthogonal unit cell basis vector lengths, in angstroms
        :multiplicity: Optional, number of unit cells in the a, b, c directions (respectively); default is (1, 1, 1)
        :overwrite: Optional, disables overwriting if False (default is False)
        :deb_dict (dict): Optional, dictionary containing the Debeye-Waller parameters for each unique element in
            the XYZ object; default is None
    """
    _overwrite_check(outfile, overwrite)
    if deb_dict is None:
        deb_dict = {}

    cel_comment = f"Number of atoms: {xyz_obj.num_sites}; {xyz_obj.comment_string}"
    # I think it would technically be possible to put everything in terms of an arbitrary metric tensor,
    #  but that seems like more effort than it's worth in almost every case
    alpha, beta, gamma = 90, 90, 90

    # Take uc_a, uc_b, uc_c in A, convert to nm
    uc_a, uc_b, uc_c = uc_params[0]/10, uc_params[1]/10, uc_params[2]/10
    a_span, b_span, c_span = (uc_a*multiplicity[0], uc_b*multiplicity[1], uc_c*multiplicity[2])
    a_min, b_min, c_min = (min([site[1] for site in xyz_obj.sites]),
                           min([site[2] for site in xyz_obj.sites]),
                           min([site[3] for site in xyz_obj.sites]))

    header = f"{'0' : <4}{a_span : ^9.6f}{b_span : ^9.6f}{c_span : ^9.6f}{alpha : ^4}{beta : ^4}{gamma : >4}\n"

    with open(outfile, 'w') as file:
        file.write(cel_comment)  # The first line in a cel file is a comment string
        file.write(header)  # The second line in a cel file defines its extent

        for site in xyz_obj.sites:
            e = site[0]
            x = f"{round((float(site[1])/10 + abs(a_min)) / a_span, 6) : .6f}"
            y = f"{round((float(site[2])/10 + abs(b_min)) / b_span, 6) : .6f}"
            z = f"{round((float(site[3])/10 + abs(c_min)) / c_span, 6) : .6f}"
            b = f"{deb_dict.get(site[0], 0.0) : .6f}"  # Default b_iso is 0

            line = f"{e : <3}{x : ^9}{y : ^9}{z : ^9}{'1.0' : ^5}{b : ^9}{' 0.0 0.0 0.0' : >11}\n"
            file.write(line)

        # EoF indicator
        file.write(r"*")


#################
# CMD Mode Only #
#################

def _parse_args(args: typing.Sequence[str]) -> Namespace:
    parser = argparse.ArgumentParser(
        description="Convert from the .xyz file format.")
    parser.add_argument("infile", nargs="+", default="",
                        help="The path to the input file; currently only supports xyz files as " +
                             "input.  Multiple files can be passed at once for batch conversion.")
    parser.add_argument("--outfile", nargs="+",
                        help="The path to the output file.  If not specified, uses the same " +
                             "path as the input file with a new extension appropriate to " +
                             "--type.  If multuple infiles are specified, do not use this argument.")
    # Add additional filetypes to choices as they're implemented
    parser.add_argument("--type", nargs="+", choices=["cel", "csv"],
                        help="The type of file to output; multiple types may be listed at once, " +
                             "or type can be inferred from file extension of --outfile if provided.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Allow file overwriting on output if the output file already " +
                             "exists (default behavior is to not allow overwriting).")
    parser.add_argument("-B", "--Biso",
                        help="The path to a file containing isometric B parameters, in the format: "
                             "<E> <B> where <E> is the (unique) element symbol and <B> is the B "
                             "parameter in A^2. Only used for .cel output format.")
    parser.add_argument("--unitcell",
                        help="Three space-separated floating-point values representing the unit cell "
                             "parameters a, b, c (with an orthogonal basis) in A. Only used for "
                             ".cel output format.")
    parser.add_argument("--multiplicity",
                        help="Three space-separated integers representing the number of unit cells "
                             "in the a, b, and c directions, repsectively. Only used for .cel output format.")
    opts = parser.parse_args(args)
    return opts


def _overwrite_check(file: str | os.PathLike | typing.TextIO, flag: bool):
    if os.path.isfile(file) and not flag:
        sys.exit(f"File {file} already exists. Use the --overwrite option if you want to overwrite it.")


def _main(opts):
    ################
    # Read Data In #
    ################
    # Handle underspecification error to do with multiple input files
    if len(opts.infile) > 1 and opts.outfile is not None:
        sys.exit("When processing multiple input files, you must not specify any output files " +
                 "as this leads to ambiguous behavior (their names will be inferred from the " +
                 "list of input files).")

    for infile in opts.infile:
        infile_name, infile_ext = os.path.splitext(infile)
        if infile_ext == "":
            infile_ext = ".xyz"
        if not infile_ext == ".xyz":
            sys.exit("This program currently only supports reading from .xyz files.")

        infile = infile_name + infile_ext
        if not os.path.isfile(infile):
            sys.exit(f"File {infile} does not exist.")

        try:
            with open(infile) as file:
                print(f"Reading data from {infile}")
                xyz_obj = read_xyz(file)
        except Exception as ex:
            print(f"Something went wrong while reading from {infile}!")
            raise ex

        ##################
        # Write Data Out #
        ##################
        if opts.outfile is None:  # Outfile(s) not specified
            if opts.type is None:
                raise RuntimeError("You must specify either output type (output file names will be inferred "
                                   "from input files), or else specify an output file name (output type will "
                                   "be inferred from the output file extension).")
            types = opts.type
            if "cel" in types:
                # Ensure that unit cell params and multiplicty have been passed if output type is .cel
                if opts.unitcell is None or opts.multiplicity is None:
                    raise RuntimeError(".cel format output requires both --unitcell and --multiplicity to be specified")
            for out_type in types:
                # With no outfile provided, infer name from infile and extension from type
                outfile = f"{infile_name}.{out_type}"
                # Call the write function associated with out_type
                call = f"write_{out_type}"

                print(f"Writing data to {outfile}")
                if out_type is "cel":
                    bdict = None
                    if opts["Biso"] is not None:
                        bdict = {}
                        with open(opts["Biso"]) as bfile:
                            lines = bfile.readlines()
                        for line in lines:
                            cols = line.split()
                            bdict[cols[0]] = cols[1]
                    globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite,
                                    uc_params=opts["unitcell"], multiplicity=opts["multiplicity"], deb_dict=bdict)
                else:
                    globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite)

        else:  # Outfile(s) specified
            for outfile in opts.outfile:
                outfile_name, outfile_ext = os.path.splitext(outfile)

                # Types explicitly specified
                if opts.type is not None:
                    types = opts.type
                    for out_type in types:
                        # If --type specified, this takes precedence over given file extensions
                        outfile = f"{outfile_name}.{out_type}"
                        call = f"write_{out_type}"
                        print(f"Writing data to {outfile}")
                        if out_type is "cel":
                            bdict = None
                            if opts["Biso"] is not None:
                                bdict = {}
                                with open(opts["Biso"]) as bfile:
                                    lines = bfile.readlines()
                                for line in lines:
                                    cols = line.split()
                                    bdict[cols[0]] = cols[1]
                            globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite, deb_dict=bdict,
                                            uc_params=opts["unitcell"], multiplicity=opts["multiplicity"])
                        else:
                            globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite)

                else:  # Types not explicitly specified (attempt to infer from outfile ext)
                    out_type = outfile_ext.replace(".", "")
                    # Error if type cannot be inferred from extension
                    if out_type not in ["xyz", "csv", "cel"] and opts.type is None:
                        raise RuntimeError("Output type could not be inferred from file extension. Either specify "
                                           "output format using --type or ensure that path(s) passed to --outfile "
                                           "end with a valid extension.")
                    call = f"write_{out_type}"
                    print(f"Writing data to {outfile}")
                    if out_type is "cel":
                        bdict = None
                        if opts["Biso"] is not None:
                            bdict = {}
                            with open(opts["Biso"]) as bfile:
                                lines = bfile.readlines()
                            for line in lines:
                                cols = line.split()
                                bdict[cols[0]] = cols[1]
                        globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite,
                                        uc_params=opts["unitcell"], multiplicity=opts["multiplicity"], deb_dict=bdict)
                    else:
                        globals()[call](xyz_obj=xyz_obj, outfile=outfile, overwrite=opts.overwrite)
        print()  # Line breaks after output loops look nice
    # Everything was successful!
    print("Done!")
    sys.exit(0)


# Run main loop iff called from the command line, but not when imported as module
if __name__ == "__main__":
    import sys
    import argparse
    opts = _parse_args(sys.argv)
    _main(opts)

#%%
