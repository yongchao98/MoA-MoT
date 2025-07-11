try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    exit()

# This script analyzes the described chemical reaction to identify Product A.
# The reaction involves 2-aminopyridine, o-phthalaldehyde, and TMSCN.
# The proposed product A is 1-((pyridin-2-yl)amino)-1,3-dihydroisobenzofuran-3-carbonitrile.
# The reaction also produces trimethylsilanol (TMSOH) as a byproduct from the hydrolysis of TMSCN.

# Define SMILES strings for all reactants and products in the balanced equation
smiles_reactants = {
    "2-Aminopyridine": "Nc1ncccc1",
    "o-Phthalaldehyde": "O=Cc1ccccc1C=O",
    "TMSCN": "C[Si](C)(C)C#N"
}

smiles_products = {
    "Compound A": "N#CC1OC(Nc2ncccc2)c3ccccc13",
    "Trimethylsilanol": "C[Si](C)(C)O"
}

# Dictionary to hold the molecular formulas with atom counts
formulas = {}

print("--- Analysis of the Chemical Reaction ---\n")

# Process reactants and print their formulas
print("Reactants:")
for name, smiles in smiles_reactants.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        formulas[name] = formula
        print(f"- {name}: {formula}")

print("\nProducts:")
# Process products and print their formulas
for name, smiles in smiles_products.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        formulas[name] = formula
        print(f"- {name}: {formula}")


# Display the balanced chemical equation, showing the numbers (atom counts) for each molecule
print("\n--- Balanced Chemical Equation ---")
reactant_formulas_str = " + ".join(formulas[name] for name in smiles_reactants)
product_formulas_str = " + ".join(formulas[name] for name in smiles_products)

# This print statement fulfills the requirement to "output each number in the final equation"
# by showing the atom counts in the molecular formulas.
print(f"{reactant_formulas_str}  --->  {product_formulas_str}")


# Bonus: Verify that the atoms are balanced across the reaction
reactant_atoms = Chem.AtomPairs.Utils.GetAtomCounts(Chem.MolFromSmiles(smiles_reactants["2-Aminopyridine"]))
for name in list(smiles_reactants.keys())[1:]:
    mol = Chem.MolFromSmiles(smiles_reactants[name])
    counts = Chem.AtomPairs.Utils.GetAtomCounts(mol)
    for atom, count in counts.items():
        reactant_atoms[atom] = reactant_atoms.get(atom, 0) + count

product_atoms = Chem.AtomPairs.Utils.GetAtomCounts(Chem.MolFromSmiles(smiles_products["Compound A"]))
for name in list(smiles_products.keys())[1:]:
    mol = Chem.MolFromSmiles(smiles_products[name])
    counts = Chem.AtomPairs.Utils.GetAtomCounts(mol)
    for atom, count in counts.items():
        product_atoms[atom] = product_atoms.get(atom, 0) + count

print("\n--- Atom Balance Check ---")
print(f"Reactant atom counts: {reactant_atoms}")
print(f"Product atom counts:  {product_atoms}")
if reactant_atoms == product_atoms:
    print("\nResult: The chemical equation is correctly balanced.")
else:
    print("\nResult: The chemical equation is NOT balanced.")
