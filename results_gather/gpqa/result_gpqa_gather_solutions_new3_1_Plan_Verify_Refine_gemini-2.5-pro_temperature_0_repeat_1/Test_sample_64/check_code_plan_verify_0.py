from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        # Add hydrogens to the molecule graph
        mol = Chem.AddHs(mol)
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        return f"Error: {e}"

def check_answer():
    """
    Checks the correctness of the proposed answer by verifying:
    1. Reaction 1: The product A must be an isomer of the reactant (conservation of atoms).
    2. Reaction 2: The product B must be a salt, not a free acid, due to the basic conditions.
    """
    # --- Data Definition ---
    # SMILES strings for all relevant compounds
    smiles_map = {
        # Reactants
        "1-vinylspiro[3.5]non-5-en-1-ol": "C=CC1(O)C2CCC=CC2C1",
        
        # Product A options
        "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "O=C1CC2=C(CCCC2)CCC1",
        "decahydro-7H-benzo[7]annulen-7-one": "O=C1C2CCCCC2CCCC1",
        
        # Product B options
        "3-ethylpent-4-enoic acid": "C=CC(CC)CC(=O)O",
        "lithium 3-ethylpent-4-enoate": "[Li+].C=CC(CC)CC(=O)[O-]"
    }

    # The proposed final answer from the LLM
    final_answer_option = 'C'
    
    # The products corresponding to the final answer
    proposed_product_A_name = "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one"
    proposed_product_B_name = "lithium 3-ethylpent-4-enoate"

    # --- Verification for Reaction 1 ---
    # A rearrangement reaction must conserve the molecular formula (i.e., the product is an isomer).
    reactant_A_formula = get_molecular_formula(smiles_map["1-vinylspiro[3.5]non-5-en-1-ol"])
    product_A_formula = get_molecular_formula(smiles_map[proposed_product_A_name])
    
    if reactant_A_formula != product_A_formula:
        # Let's also check the incorrect option to be sure
        incorrect_product_A_formula = get_molecular_formula(smiles_map["decahydro-7H-benzo[7]annulen-7-one"])
        return (f"Incorrect. The product A is not an isomer of the reactant. "
                f"Reactant formula: {reactant_A_formula}. "
                f"Proposed product A formula: {product_A_formula}. "
                f"Note: The other option for A, 'decahydro-7H-benzo[7]annulen-7-one', has formula {incorrect_product_A_formula}, which is also incorrect.")

    # --- Verification for Reaction 2 ---
    # The reaction uses a strong base (LDA) and no acidic workup is specified.
    # Therefore, the carboxylic acid product should be in its deprotonated (salt) form.
    if proposed_product_B_name == "3-ethylpent-4-enoic acid":
        return ("Incorrect. The product B is given as a free acid. "
                "Given the strong base (LDA) and lack of an acidic workup step, "
                "the product should be the deprotonated carboxylate salt (lithium 3-ethylpent-4-enoate).")
    
    # Check if the proposed product is indeed the salt
    if proposed_product_B_name != "lithium 3-ethylpent-4-enoate":
        return f"Incorrect. The name for product B, '{proposed_product_B_name}', is not the expected salt."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)
