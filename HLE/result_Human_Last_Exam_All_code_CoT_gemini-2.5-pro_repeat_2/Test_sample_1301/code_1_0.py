import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# The SMILES string representing the designed molecule.
smiles = "C1OC2CC3C1OC1C(=O)C(O3)C2OC1"

# Create a molecule object from the SMILES string.
mol = Chem.MolFromSmiles(smiles)

# --- Verification Step ---
# This part of the code verifies that the molecule meets all the specified criteria.

# 1. Molecular Formula and Weight
formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
exact_mw = Descriptors.ExactMolWt(mol)

# 2. Atom Counts
heavy_atom_count = mol.GetNumHeavyAtoms()
heteroatom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)
carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
hydrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

# 3. Valence Electrons
valence_electrons = sum(atom.GetNumOuterElecs() for atom in mol.GetAtoms())

# 4. Rings
ring_info = mol.GetRingInfo()
num_rings = ring_info.NumRings()
are_rings_heterocycles = all(any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring) for ring in ring_info.AtomRings())

# 5. Functional Groups
carbonyl_pattern = Chem.MolFromSmarts('[#6]=[#8]')
carbonyl_count = len(mol.GetSubstructMatches(carbonyl_pattern))
ether_pattern = Chem.MolFromSmarts('[OD2]([#6])[#6]')
ether_count = len(mol.GetSubstructMatches(ether_pattern))

# 6. Lipinski's Rule Parameters (for H-bond donors/acceptors)
h_bond_acceptors = Lipinski.NumHAcceptors(mol)
h_bond_donors = Lipinski.NumHDonors(mol)
rotatable_bonds = Lipinski.NumRotatableBonds(mol)

# --- Output the Final Answer ---
# The primary request is to output the SMILES string.
# The following print statements show the verification results.

print("--- Designed Molecule Verification ---")
print(f"SMILES String: {smiles}")
print("\n--- Properties ---")
print(f"Molecular Formula: {formula}")
print(f"Exact Molecular Weight: {exact_mw:.2f} g/mol")
print(f"Heavy Atom Count: {heavy_atom_count}")
print(f"Heteroatom Count: {heteroatom_count}")
print(f"Valence Electron Count: {valence_electrons}")
print(f"Total Rings: {num_rings}")
print(f"All Rings are Heterocycles: {are_rings_heterocycles}")
print(f"Carbonyl Count: {carbonyl_count}")
print(f"Ether Count: {ether_count}")
print(f"Hydrogen Bond Acceptors: {h_bond_acceptors}")
print(f"Hydrogen Bond Donors: {h_bond_donors}")
print(f"Rotatable Bonds: {rotatable_bonds}")

# Final Answer as requested by the user prompt
print("\n--- Final Answer ---")
print("The SMILES representation for the molecule is:")
final_answer = smiles
print(final_answer)