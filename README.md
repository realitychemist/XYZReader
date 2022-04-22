# XYZReader

***A module to work with xyz-formatted files.***

This module has callable functions to read xyz-formatted atomic position files, as well as to write them.  In addition, it can be used as a standalone script for converting filetypes.

## Usage

### Import Mode

There are two functions that can be used in import mode: `readXYZ` and `writeXYZ`, whose purpose should be pretty self-evident.

#### `readXYZ(file)`
`readXYZ` expects a properly-formatted `.xyz` file, and doesn't do any of the filesystem checking work.  You should handle things like making sure the file exists in your own code.

Since there are several fairly common variations on the `.xyz` file structure, the meaning of properly-formatted is:
1. The first line contains an integer (the number of coordinate lines in the file)
2. The second line contains some kind of comment string describing the file
3. The rest of the lines in the file contain four space-separated entries: an atomic number or symbol, and the x, y, and z coordinates of that atom
This is the style of `.xyz` file supported by the UCSF software [Chimera](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/xyz.html).

`readXYZ` will return an instance of the `XYZ` object, which has several fields:
* `num_sites` contains the number of sites described in the file
* `comment_string` is the comment string from the second line of the `.xyz` file, as a python string
* `sites` is a list of lists, where each sublist describes one site from the file as `[E, x, y, z]`
  * `E` is a string containing the atomic number or symbol associated with the site
  * `x`, `y`, and `z` are the coordinates of that site; units depend on the input file, no unit conversion is done
  
 #### `writeXYZ(xyz_obj, outfile)`
 `writeXYZ` takes an instance of an `XYZ` object and a path to where you'd like to write the file (as a string), as described above, and writes the data contained in it to an `.xyz` file (formatted as described above).  Just like `readXYZ`, no filesystem checking is performed here.  Make sure the location you've specified is writable.
 
 ### Script Mode
 Script mode is meant to be used for in-place conversion to and from `.xyz` files.  It is currently early days for this part of the module.  Currently, it only converts *from* `.xyz`, and supports two output formats: `.csv` and `.cel`.
 
 ## Planned Features
 * Output to more file formats (`.cif`, maybe z-matrix)
 * Convert other filetypes to `.xyz`
