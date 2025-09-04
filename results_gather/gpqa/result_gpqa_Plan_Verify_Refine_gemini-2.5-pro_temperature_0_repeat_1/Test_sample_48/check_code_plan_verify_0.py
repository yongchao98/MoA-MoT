import sys
import re

# Try to import the RDKit library, which is essential for this chemical analysis.
# If not found, provide instructions to install it.
try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    # If rdkit is not installed, the check cannot be performed.
    # We print an informative error message instead of the final result.
    # This is a common practice for scripts with external dependencies.
    print("Error: The 'rdkit' library is required to run this check. Please install it using 'pip install rdkit-pypi'.")
    sys.exit(1) # Exit the script with an error code.

def get_formula_from_smiles(smiles: str) -> str:
    """
    Converts a SMILES string into its molecular formula using RDKit.
    Returns None if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return CalcMolFormula(mol)

def parse_formula(formula_str: str) -> dict:
    """
    Parses a molecular formula string (e.g., 'C8H10') into a dictionary of atom counts.
    """
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    atoms = {}
    for element, count in pattern.findall(formula_str):
        atoms[element] = atoms.get(element, 0) + (int(count) if count else 1)
    return atoms

def subtract_formulas(formula_total: str, formula_product: str) -> str:
    """
    Calculates the difference between two molecular formulas.
    Used to find the formula of the eliminated molecule(s) in a reaction.
    """
    atoms_total = parse_formula(formula_total)
    atoms_product = parse_formula(formula_product)
    diff = {}
    for atom, count in atoms_total.items():
        diff[atom] = count - atoms_product.get(atom, 0)
    
    # Reconstruct the formula string for the difference
    diff_str = "".join([f"{atom}{count if count > 1 else ''}" for atom, count in sorted(diff.items()) if count > 0])
    return diff_str

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the chemical reactions' stoichiometry.
    """
    llm_selected_option = 'B'

    # Define reactants and products using SMILES for unambiguous chemical representation.
    # These SMILES are derived from the IUPAC names in the question.
    reactants = {
        'A1': 'CC(OC)(OC)N',  # 1,1-dimethoxyethan-1-amine
        'A2': 'C=CC(C)O',     # but-3-en-2-ol
        'B': 'C#C[C@H](C)[C@@H](C)C#C', # (3R,4S)-3,4-dimethylhexa-1,5-diyne
        'C': 'C=C(CC)COOC=C' # 2-((vinyloxy)methyl)but-1-ene
    }

    options = {
        'A': {
            'A': 'N/C=C/OC(=C(C))C',
            'B': 'C/C=C1/C=C/C1=C\\C',
            'C': 'OCCCC(=C)CC'
        },
        'B': {
            'A': 'NC1CCOC=C1C',      # 6-methyl-3,4-dihydro-2H-pyran-2-amine
            'B': 'C/C=C1\\CCC/C1=C/C', # (1Z,2E)-1,2-diethylidenecyclobutane
            'C': 'O=CCCC(=C)CC'      # 4-methylenehexanal
        },
        'C': {
            'A': 'N/C=C/OC(=C(C))C',
            'B': 'C/C=C1/C=C/C1=C\\C',
            'C': 'O=CCCC(=C)CC'
        },
        'D': {
            'A': 'NC1CCOC=C1C',
            'B': 'C/C=C1\\CCC/C1=C/C',
            'C': 'OCCCC(=C)CC'
        }
    }
    
    # --- Verification Step 1: Check Reaction C (Claisen Rearrangement) ---
    # This is an intramolecular rearrangement, so mass must be conserved.
    # Product formula must equal reactant formula.
    reactant_c_formula = get_formula_from_smiles(reactants['C'])
    product_c_formula = get_formula_from_smiles(options[llm_selected_option]['C'])
    
    if reactant_c_formula != product_c_formula:
        return (f"Incorrect. Constraint for Reaction C is not satisfied. "
                f"This is an intramolecular rearrangement, so formulas should match. "
                f"Reactant C formula: {reactant_c_formula}, "
                f"Proposed Product C formula: {product_c_formula}.")

    # --- Verification Step 2: Check Reaction B (Cope Rearrangement) ---
    # This is also an intramolecular rearrangement. Mass must be conserved.
    reactant_b_formula = get_formula_from_smiles(reactants['B'])
    product_b_formula = get_formula_from_smiles(options[llm_selected_option]['B'])

    if reactant_b_formula != product_b_formula:
        return (f"Incorrect. Constraint for Reaction B is not satisfied. "
                f"This is an intramolecular rearrangement, so formulas should match. "
                f"Reactant B formula: {reactant_b_formula}, "
                f"Proposed Product B formula: {product_b_formula}.")

    # As a sanity check, verify the LLM's reasoning that the alternative product for B is wrong.
    product_b_alt_formula = get_formula_from_smiles(options['C']['B'])
    if product_b_alt_formula == reactant_b_formula:
        return (f"Incorrect. The LLM's reasoning is flawed. It claimed the product B from option C "
                f"({product_b_alt_formula}) has the wrong formula, but it matches the reactant ({reactant_b_formula}).")

    # --- Verification Step 3: Check Reaction A (Condensation) ---
    # This reaction involves eliminating small molecules (H+ and Heat suggest this).
    # We check if the mass difference corresponds to a plausible elimination.
    reactant_a1_formula = get_formula_from_smiles(reactants['A1'])
    reactant_a2_formula = get_formula_from_smiles(reactants['A2'])
    product_a_formula = get_formula_from_smiles(options[llm_selected_option]['A'])
    
    # Combine formulas of the two reactants
    total_reactant_a_atoms = parse_formula(reactant_a1_formula)
    reactant_a2_atoms = parse_formula(reactant_a2_formula)
    for atom, count in reactant_a2_atoms.items():
        total_reactant_a_atoms[atom] = total_reactant_a_atoms.get(atom, 0) + count
    total_reactant_a_formula = "".join([f"{atom}{count if count > 1 else ''}" for atom, count in sorted(total_reactant_a_atoms.items())])

    # The reaction eliminates two molecules of methanol (CH3OH), so the formula of the eliminated part should be C2H8O2.
    expected_elimination = "C2H8O2"
    actual_elimination = subtract_formulas(total_reactant_a_formula, product_a_formula)

    if actual_elimination != expected_elimination:
        return (f"Incorrect. Constraint for Reaction A is not satisfied. "
                f"The stoichiometry is incorrect. The difference between reactants ({total_reactant_a_formula}) "
                f"and product ({product_a_formula}) is {actual_elimination}, "
                f"which does not match the expected elimination of two methanol molecules ({expected_elimination}).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the final result.
result = check_answer_correctness()
print(result)