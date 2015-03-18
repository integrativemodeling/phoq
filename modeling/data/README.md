`angle.py` and `dihedral.py` generate the data files for the CA forcefield
used in this modeling (`CAAngleRestraint.dat` and `CADihedralRestraint.dat`
respectively). They take as input a file `file.list` which is a list of
filenames, each of which is the output of DSSP of a protein (typically
`file.list` will cover a large subset of PDB).

The format of the angle data file is simple; the first column lists the
secondary structure class of a set of three contiguous residues, the second
the angle between the CA atoms in those three residues (binned into 1 degree
bins), and the third the number of times this is seen in PDB. The dihedral
data file is similar, except that it lists five residues, the dihedral
angle between the first four CA atoms, the dihedral between the last four
(both in 10 degree bins), and the frequency.

These data files are parsed by the modeling scripts (there is also similar, but
not quite identical, code in the
IMP.pmi.restraints.stereochemistry.SecondaryStructure class) and used to
create IMP.atom.CAAngleRestraint and IMP.atom.CADihedralRestraint objects.
