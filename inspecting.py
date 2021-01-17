from MDAnalysis import Universe


obj = Universe("open_V_ps.pdb","open_450ns_protein_only.dcd")
obj_ref = Universe("MUTATED.pdb")
obj_ref2 = Universe("open_V_ps.pdb")

print(len(obj_ref.residues))
print(len(obj_ref2.residues))

"""
for res in obj_ref.atoms.posistions:
    positions = obj_ref.atoms.positions
    print(positions)"""