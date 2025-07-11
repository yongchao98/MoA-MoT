from rdkit import Chem
from rdkit.Chem import Descriptors

# The user's request involves identifying the heavier product from a chemical reaction.
# The starting material's SMILES string 'CC12COC(OC1)(OC2)C1=CC=CC=C1' is malformed,
# so we must deduce the reactants based on chemical principles and clues in the string.
# The plan is to identify the two hydrolysis products (a ketone and a diol) based on these clues,
# calculate their molar masses, and return the SMILES string of the heavier one.

# Product 1: Acetophenone, deduced from 'CC' and 'C1=CC=CC=C1' fragments.
ketone_smiles = "CC(=O)c1ccccc1"
ketone_mol = Chem.MolFromSmiles(ketone_smiles)
ketone_mw = Descriptors.MolWt(ketone_mol)

# Product 2: Diethylene glycol, deduced from the presence of 3 oxygen atoms and a 'COC' ether link.
diol_smiles = "OCCOCCO"
diol_mol = Chem.MolFromSmiles(diol_smiles)
diol_mw = Descriptors.MolWt(diol_mol)

print("The acid-catalyzed hydrolysis reaction produces two main products:")
print(f"1. Ketone (Acetophenone), SMILES: {ketone_smiles}")
print(f"   Molar Mass of Acetophenone = {ketone_mw:.2f} g/mol")
print(f"2. Diol (Diethylene Glycol), SMILES: {diol_smiles}")
print(f"   Molar Mass of Diethylene Glycol = {diol_mw:.2f} g/mol")

print("\nComparing the two molar masses:")
if ketone_mw > diol_mw:
    print(f"Acetophenone ({ketone_mw:.2f}) > Diethylene Glycol ({diol_mw:.2f})")
    print("Therefore, the product with the higher molar mass is Acetophenone.")
    final_answer_smiles = ketone_smiles
else:
    print(f"Diethylene Glycol ({diol_mw:.2f}) > Acetophenone ({ketone_mw:.2f})")
    print("Therefore, the product with the higher molar mass is Diethylene Glycol.")
    final_answer_smiles = diol_smiles

print("\nThe SMILES string for the product with the higher molar mass is:")
print(final_answer_smiles)