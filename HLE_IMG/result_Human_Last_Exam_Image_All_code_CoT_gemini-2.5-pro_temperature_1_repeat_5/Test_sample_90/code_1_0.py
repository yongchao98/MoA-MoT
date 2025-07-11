from rdkit import Chem
from rdkit.Chem import Draw

# Define the SMILES strings for the products A, B, and C
smiles_a = "CC(=O)N1C2C=C(C(=O)OC)CN2C1(CCCC1)C1=NCCC1" # Product A
smiles_b = "COC(=O)/C=C/C1N(C(=O)C2=NCCC2)CCCC1C(=O)O" # A plausible structure for Product B
smiles_c = "CC(=O)N1CCCC1(C(=O)O)C2=NCCC2"          # Product C

# Generate molecule objects from SMILES
mol_a = Chem.MolFromSmiles(smiles_a)
mol_b = Chem.MolFromSmiles(smiles_b)
mol_c = Chem.MolFromSmiles(smiles_c)

# Create a list of molecules and labels for drawing
mols = [mol_c, mol_a, mol_b]
labels = ["Product C", "Product A", "Product B"]

# Add labels to molecules
for mol, label in zip(mols, labels):
    mol.SetProp("_Name", label)

# Draw the molecules in a grid
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300, 300), legends=[mol.GetProp("_Name") for mol in mols])

# To display the image in environments like Jupyter notebooks
# In a standard python script, we save it to a file.
img.save("reaction_products.png")

print("The structures of the three products are:")
print(f"Product A: {Chem.MolToSmiles(mol_a)}")
print(f"Product B: {Chem.MolToSmiles(mol_b)}")
print(f"Product C: {Chem.MolToSmiles(mol_c)}")
print("\nAn image of the structures has been saved as 'reaction_products.png'")