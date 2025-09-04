# First, ensure you have rdkit installed:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def check_correctness():
    """
    Checks the correctness of the answer by verifying the chemical logic.
    1. Checks the molecular formula of all options.
    2. Checks for the presence of the key "p-tolyl" substructure.
    3. Verifies that the chosen answer is the unique option satisfying the constraints.
    """
    # --- Data Definition ---
    # The options as defined in the question
    options = {
        'A': {'name': '2-(1-phenylprop-1-en-2-yl)oxirane', 'smiles': 'CC(=Cc1ccccc1)C1CO1'},
        'B': {'name': '2-methyl-3-styryloxirane', 'smiles': 'CC1OC(C=Cc2ccccc2)C1'},
        'C': {'name': '2-styrylepoxide', 'smiles': 'c1ccc(cc1)C=CC1CO1'},
        'D': {'name': '2-(4-methylstyryl)oxirane', 'smiles': 'Cc1ccc(C=CC2CO2)cc1'}
    }

    # The product identified from NMR data in the reasoning
    product = {'name': '4-(4-methylphenyl)but-3-en-2-one', 'smiles': 'CC(=O)C=Cc1ccc(C)cc1'}

    # The final answer to be checked
    final_answer = 'D'

    # --- Constraints from the Question and Reasoning ---
    required_formula = "C11H12O"
    # SMARTS pattern for a p-tolyl group (a methyl-substituted benzene ring with another attachment)
    p_tolyl_smarts = '[#6]1(-[#6])ccc(-[#6&H3])cc1'

    # --- Helper Functions ---
    def get_formula(smiles_str):
        """Calculates the molecular formula from a SMILES string."""
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return "Invalid SMILES"
        return CalcMolFormula(mol)

    def has_substructure(smiles_str, smarts_pattern):
        """Checks if a molecule contains a given substructure."""
        mol = Chem.MolFromSmiles(smiles_str)
        pattern = Chem.MolFromSmarts(smarts_pattern)
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)

    # --- Verification Logic ---

    # 1. Verify the properties of the identified product (validates the reasoning's starting point)
    product_formula = get_formula(product['smiles'])
    if product_formula != required_formula:
        return f"Reasoning is flawed: The identified product '{product['name']}' has formula {product_formula}, not {required_formula}."
    if not has_substructure(product['smiles'], p_tolyl_smarts):
        return f"Reasoning is flawed: The identified product '{product['name']}' does not contain a p-tolyl group."

    # 2. Find all candidate options that meet the necessary criteria
    valid_candidates = []
    for key, data in options.items():
        formula = get_formula(data['smiles'])
        has_tolyl = has_substructure(data['smiles'], p_tolyl_smarts)
        
        # A valid candidate must have the correct formula AND the p-tolyl group
        if formula == required_formula and has_tolyl:
            valid_candidates.append(key)

    # 3. Check the final answer against the valid candidates
    if final_answer not in valid_candidates:
        chosen_option_data = options[final_answer]
        chosen_formula = get_formula(chosen_option_data['smiles'])
        if chosen_formula != required_formula:
            return f"Incorrect. The chosen answer '{final_answer}' is wrong because its molecular formula is {chosen_formula}, but the question requires {required_formula}."
        if not has_substructure(chosen_option_data['smiles'], p_tolyl_smarts):
            return f"Incorrect. The chosen answer '{final_answer}' is wrong because it does not contain the required p-tolyl group, which is present in the product."
        # This part should not be reached if logic is sound, but as a fallback:
        return f"Incorrect. The chosen answer '{final_answer}' does not meet the required criteria."

    if len(valid_candidates) > 1:
        return f"Incorrect. The reasoning that only one option is valid is flawed. Options {valid_candidates} all have the correct formula and contain a p-tolyl group."

    # 4. If the final answer is the *only* valid candidate, the logic is sound.
    if valid_candidates == [final_answer]:
        return "Correct"
    
    # Fallback for any unhandled case
    return "Could not determine correctness."

# Run the check
result = check_correctness()
print(result)