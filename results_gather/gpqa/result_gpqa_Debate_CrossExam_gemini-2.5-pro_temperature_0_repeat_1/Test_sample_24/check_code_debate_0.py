import sys
from collections import Counter
import re

# As rdkit is a specialized library, we handle its potential absence.
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit library not found. Cannot perform the chemical validation.")
    print("Please install it using: pip install rdkit")
    # Exit gracefully if the required library is not available.
    sys.exit(0)

def get_atom_counts(formula_str: str) -> Counter:
    """
    Parses a molecular formula string (e.g., 'C12H22O2') into a Counter
    object mapping each element to its count.
    """
    atom_counts = Counter()
    # Regex to find element symbols (e.g., C, H, O, Cl) and their optional counts
    for element, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula_str):
        atom_counts[element] += int(count) if count else 1
    return atom_counts

def check_answer():
    """
    Checks the correctness of the proposed answer by verifying the
    stoichiometry of the two chemical reactions.
    """
    # SMILES strings for the reactants from the proposed answer (Option A) and the products from the question.
    # Note: A specific isomer is chosen for reactant A that fits the name. The molecular formula is the key check.
    smiles_map = {
        "Reactant A": "CC1CCC2CCCC(C)C2C1(O)O",  # An isomer of 2,7-dimethyloctahydronaphthalene-4a,8a-diol
        "Product 1": "CC1CCC2(CC(C)CC2)C(=O)C1", # 2,8-dimethylspiro[4.5]decan-6-one
        "Reactant B": "c1ccccc1COCC=C(C)C",      # (((3-methylbut-2-en-1-yl)oxy)methyl)benzene
        "Product 2": "CC(=CCC(c1ccccc1)O)C"      # 4-methyl-1-phenylpent-3-en-1-ol
    }

    formulas = {}
    for name, smiles in smiles_map.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return f"Error: Could not parse SMILES string for '{name}'."
        formulas[name] = get_atom_counts(CalcMolFormula(mol))

    # --- Check Reaction 1: Pinacol Rearrangement (A -> Product 1 + H2O) ---
    # The reaction involves the loss of one water molecule.
    # So, atoms in A should equal atoms in Product 1 plus atoms in H2O.
    expected_reactant_A_atoms = formulas["Product 1"] + Counter({'H': 2, 'O': 1})
    if formulas["Reactant A"] != expected_reactant_A_atoms:
        return (f"Incorrect: Stoichiometry for Reaction 1 (Pinacol rearrangement) is wrong.\n"
                f"Reactant A formula: {formulas['Reactant A']}\n"
                f"Expected formula (Product 1 + H2O): {expected_reactant_A_atoms}")

    # --- Check Reaction 2: [2,3]-Wittig Rearrangement (B -> Product 2) ---
    # This is an intramolecular rearrangement, so atom counts must be identical.
    if formulas["Reactant B"] != formulas["Product 2"]:
        return (f"Incorrect: Stoichiometry for Reaction 2 ([2,3]-Wittig rearrangement) is wrong.\n"
                f"Reactant B formula: {formulas['Reactant B']}\n"
                f"Product 2 formula: {formulas['Product 2']}\n"
                f"The formulas should be identical for a rearrangement.")

    # If both checks pass, the answer is correct from a stoichiometric and reaction-type standpoint.
    return "Correct"

# Run the check and print the result.
result = check_answer()
print(result)