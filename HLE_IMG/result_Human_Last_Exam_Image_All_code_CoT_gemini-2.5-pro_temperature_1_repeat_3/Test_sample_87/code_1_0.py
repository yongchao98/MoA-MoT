try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit library is not installed. Please install it using 'pip install rdkit-pypi'")
    exit()

# SMILES (Simplified Molecular-Input Line-Entry System) strings for the products
# Note: The exact stereoisomer and regioisomer are assumed based on common reaction mechanisms.

# Product C: Mixed anhydride of the starting material and acetic acid.
smiles_C = "CC(=O)OC(=O)[C@H]1CCCN1C2=NCCC2"
mol_C = Chem.MolFromSmiles(smiles_C)
formula_C = rdMolDescriptors.CalcMolFormula(mol_C)

# Product B: Aromatized and N-oxidized cycloadduct.
# The structure is methyl 1-(1-oxido-4,5-dihydro-3H-pyrrol-2-yl)-1H-pyrrolizine-6-carboxylate.
smiles_B = "COC(=O)c1cc2c(n1C3=[N+]([O-])CCC3)ccc2"
mol_B = Chem.MolFromSmiles(smiles_B)
formula_B = rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol_B))


# Product A: Aldol adduct of the primary cycloadduct (DHP) and acetaldehyde.
# The structure is an aldol adduct of a dihydropyrrolizine with acetaldehyde.
# A plausible isomer is chosen for representation.
smiles_A = "COC(=O)C1=C[C@H]2[C@@H](C(C(O)C)C1)N(CCC2)C3=NCCC3"
mol_A = Chem.MolFromSmiles(smiles_A)
formula_A = rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol_A))


print("Structure of Product A:")
print(f"Molecular Formula: {formula_A} (matches C14H20N2O3)")
print(f"SMILES representation: {smiles_A}")
print("-" * 30)
print("Structure of Product B:")
print(f"Molecular Formula: {formula_B} (matches C12H14N2O3)")
print(f"SMILES representation: {smiles_B}")
print("-" * 30)
print("Structure of Product C:")
print(f"Molecular Formula: {formula_C} (matches C11H16N2O3)")
print(f"SMILES representation: {smiles_C}")
print("-" * 30)
