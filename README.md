Convert between fractional crystal coordinates and Cartesian crystal coordinates.
Supports input file of VASP \(POSCAR\) and input file of pw.x of Quantum ESPRESSO.

Usage: Cart2Dire.exe  \[QE_pw_input_file\]

If QE_pw_input_file is not provided, 
it will try to read the POSCAR of VASP in the current directory, 
and if POSCAR is not found, it will ask for the name of the Quantum ESPRESSO 
pw.x input file.
The output name has the same suffix \(if any\) with the input file, and the file name 
is appended by \_new.
if the input file is fractional crystal coordinates \(a.k.a. direct coordinates\), 
the output file will be Cartesian coordinates, and vice versa.

