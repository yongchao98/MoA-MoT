import sys
from collections import Counter

# Try to import RDKit. If it's not available, provide a helpful message.
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    # If RDKit is not installed, the check cannot be performed.
    # We will print an error message and exit.
    print("Error: RDKit library not found.")
    print("Please install it to run this check: pip install rdkit")
    sys.exit("RDKit is a required dependency for this verification script.")

def get_formula_from_smiles(smiles_string):
    """Converts a SMILES string to a molecular formula."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return f"Invalid SMILES: {smiles_string}"
    return rdMolDescriptors.CalcMolFormula(mol)

def get_combined_formula_from_smiles(smiles_list):
    """Calculates a combined molecular formula from a list of SMILES strings."""
    total_counts = Counter()
    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return f"Invalid SMILES in list: {s}"
        # Use RDKit's formula calculation and parse it
        formula = rdMolDescriptors.CalcMolFormula(mol)
        import re
        atoms = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for atom, count in atoms:
            total_counts[atom] += int(count) if count else 1
    
    # Format the formula string in Hill order (C, H, then alphabetical)
    formula_str = ""
    if 'C' in total_counts:
        formula_str += f"C{total_counts.pop('C') if total_counts['C'] > 1 else ''}"
    if 'H' in total_counts:
        formula_str += f"H{total_counts.pop('H') if total_counts['H'] > 1 else ''}"
    
    for atom in sorted(total_counts.keys()):
        formula_str += f"{atom}{total_counts[atom] if total_counts[atom] > 1 else ''}"
        
    return formula_str

def check_reactions():
    """
    Checks the correctness of the products for the three reactions
    as proposed in answer B.
    """
    errors = []

    # --- Reaction C: Claisen Rearrangement ---
    # Reactant: 2-((vinyloxy)methyl)but-1-ene -> Product C: 4-methylenehexanal
    # This is a [3,3]-sigmatropic rearrangement, so the product must be an isomer.
    smiles_c_reactant = "C=C(CC)COC=C"
    smiles_c_product = "O=CCCC(=C)CC"
    formula_c_reactant = get_formula_from_smiles(smiles_c_reactant)
    formula_c_product = get_formula_from_smiles(smiles_c_product)

    if formula_c_reactant != formula_c_product:
        errors.append(
            f"Constraint Violated (Reaction C): Product C is not an isomer of its reactant.\n"
            f"  - Reactant Formula: {formula_c_reactant}\n"
            f"  - Proposed Product Formula: {formula_c_product}"
        )
    # Mechanistic check: The product of a Claisen rearrangement of an allyl vinyl ether
    # is a gamma,delta-unsaturated carbonyl. 4-methylenehexanal fits this description.

    # --- Reaction B: Hopf Rearrangement ---
    # Reactant: (3R,4S)-3,4-dimethylhexa-1,5-diyne -> Product B: (3Z,4E)-3,4-diethylidenecyclobut-1-ene
    # This is also a rearrangement, so the product must be an isomer.
    smiles_b_reactant = "C#C[C@H](C)[C@@H](C)C#C"
    smiles_b_product = "C/C=C1/C(=C\C)C=C1" # SMILES for (3Z,4E)-3,4-diethylidenecyclobut-1-ene
    formula_b_reactant = get_formula_from_smiles(smiles_b_reactant)
    formula_b_product = get_formula_from_smiles(smiles_b_product)

    if formula_b_reactant != formula_b_product:
        errors.append(
            f"Constraint Violated (Reaction B): Product B is not an isomer of its reactant.\n"
            f"  - Reactant Formula: {formula_b_reactant}\n"
            f"  - Proposed Product Formula: {formula_b_product}"
        )
    # Mechanistic check: The thermal rearrangement of a 1,5-diyne (Hopf rearrangement)
    # produces a cyclobutene derivative, not a cyclobutane. The proposed product is a
    # cyclobutene, which is correct. This check is crucial for eliminating other options.

    # --- Reaction A: Eschenmoser-Claisen type reaction ---
    # Reactants: 1,1-dimethoxyethan-1-amine + but-3-en-2-ol -> Product A + 2 MeOH
    # This is a condensation reaction. We must check the atom balance.
    smiles_a_reactants = ["N[CH](OC)OC", "C=CC(C)O"]
    smiles_a_product_and_byproducts = ["N/C=C/OC(C)=CC", "CO", "CO"] # Product A + 2 Methanol
    
    formula_a_reactants = get_combined_formula_from_smiles(smiles_a_reactants)
    formula_a_products = get_combined_formula_from_smiles(smiles_a_product_and_byproducts)

    if formula_a_reactants != formula_a_products:
        errors.append(
            f"Constraint Violated (Reaction A): Atom balance is incorrect.\n"
            f"  - Combined Reactant Formula: {formula_a_reactants}\n"
            f"  - Combined Product Formula (Product A + 2 MeOH): {formula_a_products}"
        )

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. The following constraints were not satisfied:\n" + "\n".join(errors)

# Run the check and print the result.
result = check_reactions()
print(result)