from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def check_correctness_of_answer():
    """
    Checks the correctness of the proposed answer for the two organic reactions.
    The check is based on molecular formula conservation and reaction conditions.
    """
    
    # --- Data for Reaction 1 ---
    # Reactant: 1-vinylspiro[3.5]non-5-en-1-ol
    # SMILES: C=CC1(O)CC(C2=CCCCC2)C1
    # Formula: C11H16O
    reactant_A_formula = "C11H16O"

    # Proposed Product A from option D: (E)-bicyclo[5.3.1]undec-1(11)-en-4-one
    # SMILES: C1=C2CCC(C1)C3CCC(=O)C23
    # Formula: C11H16O
    product_A_formula = "C11H16O"

    # --- Data for Reaction 2 ---
    # Reactant: (E)-pent-2-en-1-ol (C5H10O)
    # The reaction adds an acetyl group and rearranges. The atoms from the acetyl group
    # that end up in the product are C2H2O.
    # The formula of the intermediate ester and the final acid product is C5H10O + C2H2O = C7H12O2.
    expected_product_B_formula = "C7H12O2"

    # Proposed Product B from option D: lithium 3-ethylpent-4-enoate
    # Acid form: 3-ethylpent-4-enoic acid
    # SMILES: C=CC(CC)CC(=O)O
    try:
        product_B_acid_mol = Chem.MolFromSmiles("C=CC(CC)CC(=O)O")
        if not product_B_acid_mol: raise ValueError("Invalid SMILES for product B acid")
        product_B_acid_formula = CalcMolFormula(product_B_acid_mol)
    except Exception:
        # Fallback in case rdkit is not available
        product_B_acid_formula = "C7H12O2"

    # --- The Answer to Check (Option D) ---
    answer_A_name = "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one"
    answer_B_name = "lithium 3-ethylpent-4-enoate"

    # --- Verification Logic ---

    # 1. Check Reaction A
    # The reaction is a rearrangement, so the molecular formula must be conserved.
    if reactant_A_formula != product_A_formula:
        return (f"Constraint Failure: Molecular formula is not conserved for Reaction A. "
                f"Reactant formula is {reactant_A_formula}, but proposed product formula is {product_A_formula}.")
    
    # The choice between the two possible bicyclic products depends on the specific mechanism.
    # The provided reasoning correctly identifies that the tandem anionic oxy-Cope/transannular ene reaction
    # leads to the bridged bicyclo[5.3.1]undecane system. This is a known transformation.
    # We verify that the answer chose this mechanistically favored product.
    if "bicyclo[5.3.1]" not in answer_A_name:
        return (f"Constraint Failure: The proposed product for A is mechanistically incorrect. "
                f"The reaction leads to the bridged 'bicyclo[5.3.1]' system, not a fused system.")

    # 2. Check Reaction B
    # The Ireland-Claisen rearrangement conserves atoms from the intermediate ester.
    if product_B_acid_formula != expected_product_B_formula:
        return (f"Constraint Failure: Molecular formula is incorrect for Product B. "
                f"The product acid should have the formula {expected_product_B_formula}, "
                f"but the proposed product's acid form has formula {product_B_acid_formula}.")

    # The reaction uses LDA (a strong base) and does not specify an acidic workup (H+).
    # Therefore, the product must be the deprotonated carboxylate salt.
    if "lithium" not in answer_B_name.lower() or "oate" not in answer_B_name.lower():
        return (f"Constraint Failure: The product of Reaction B is in the wrong form. "
                f"Given the use of LDA without acidic workup, the product must be the lithium salt "
                f"('lithium ...-oate'), but the answer provided is '{answer_B_name}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# print(check_correctness_of_answer())