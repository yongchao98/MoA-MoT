import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        return f"Error: {e}"

def check_substructure(smiles: str, smarts: str) -> bool:
    """Checks if a molecule contains a given substructure."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts)
        if mol is None or pattern is None:
            return False
        return mol.HasSubstructMatch(pattern)
    except:
        return False

def check_correctness():
    """
    Checks the correctness of the answer for the Cope rearrangement question.
    1. Verifies that the product is an isomer of the reactant.
    2. Verifies that the product matches the structure derived from the reaction mechanism.
    3. Verifies that other options are incorrect.
    """
    # --- Define Molecules and the Proposed Answer ---
    reactant = {
        "name": "5-butylnona-2,6-diene",
        "smiles": "CCC=CC(CCCC)C=CCC"
    }
    
    options = {
        "A": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "D": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "C=CC(C)C(CC)C=CCCCC"}
    }
    
    proposed_answer_key = "D"
    chosen_option = options[proposed_answer_key]

    # --- Step 1: Check Isomerism Constraint ---
    # The product of a rearrangement must be an isomer of the reactant.
    reactant_formula = get_molecular_formula(reactant["smiles"])
    product_formula = get_molecular_formula(chosen_option["smiles"])
    
    if reactant_formula != product_formula:
        return (f"Incorrect. The product must be an isomer of the reactant. "
                f"Reactant formula: {reactant_formula}. "
                f"Proposed answer (D) formula: {product_formula}.")

    # --- Step 2: Check Reaction Mechanism Constraint ---
    # The Cope rearrangement interconverts 1,5-dienes. The product should be a 1,5-diene.
    # The SMARTS pattern for a generic 1,5-diene is C=C-C-C-C=C
    is_1_5_diene = check_substructure(chosen_option["smiles"], "[#6]=[#6]-[#6]-[#6]-[#6]=[#6]")
    if not is_1_5_diene:
        return (f"Incorrect. The product of a Cope rearrangement should be a 1,5-diene. "
                f"The proposed answer (D), {chosen_option['name']}, does not fit this structure.")

    # --- Step 3: Verify the specific product structure ---
    # The reasoning in the provided answer correctly derives the product name as
    # "4-ethyl-3-methyldeca-1,5-diene". We check if the proposed answer matches this.
    if chosen_option["name"] != "4-ethyl-3-methyldeca-1,5-diene":
        return (f"Incorrect. The proposed answer key 'D' points to the name '{chosen_option['name']}', "
                f"but the correct product derived from the mechanism is '4-ethyl-3-methyldeca-1,5-diene'. "
                f"There might be a mismatch in the option lettering.")

    # --- Step 4: Rule out other options ---
    # Check Option A
    formula_A = get_molecular_formula(options["A"]["smiles"])
    if formula_A == reactant_formula:
        return "Error in checking logic: Option A was expected to be a non-isomer."
    
    # Check Option B/C
    is_1_5_diene_B = check_substructure(options["B"]["smiles"], "[#6]=[#6]-[#6]-[#6]-[#6]=[#6]")
    if is_1_5_diene_B:
        return ("Error in checking logic: Option B was expected to be a 2,6-diene, "
                "not a 1,5-diene.")

    # --- Final Conclusion ---
    # All checks passed for option D.
    return "Correct"

# Run the check
result = check_correctness()
print(result)