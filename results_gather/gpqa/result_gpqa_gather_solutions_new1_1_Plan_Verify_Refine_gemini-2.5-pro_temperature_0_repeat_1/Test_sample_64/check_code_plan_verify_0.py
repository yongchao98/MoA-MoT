import re

# RDKit is a standard cheminformatics library. We'll try to use it for accuracy.
# If it's not installed, we'll use pre-calculated values as a fallback.
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not found. Using pre-calculated molecular formulas. For a more robust check, please install rdkit (`pip install rdkit`).")

def get_molecular_formula(name, smiles_map):
    """
    Calculates the molecular formula for a given chemical name.
    Uses RDKit if available, otherwise falls back to a pre-calculated dictionary.
    """
    if RDKIT_AVAILABLE:
        smiles = smiles_map.get(name)
        if not smiles:
            return f"Unknown SMILES for {name}"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return f"Invalid SMILES for {name}"
        return rdMolDescriptors.CalcMolFormula(mol)
    else:
        # Fallback to pre-calculated values
        formulas = {
            "1-vinylspiro[3.5]non-5-en-1-ol": "C11H16O",
            "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "C11H16O",
            "decahydro-7H-benzo[7]annulen-7-one": "C11H18O",
            "3-ethylpent-4-enoic acid": "C7H12O2",
            "lithium 3-ethylpent-4-enoate": "C7H11LiO2"
        }
        return formulas.get(name, f"Unknown formula for {name}")

def check_answer_correctness():
    """
    Checks the correctness of the provided answer by verifying the chemical logic.
    """
    # --- Define the problem and the proposed answer ---
    question = {
        "reaction1_reagents": "(THF, KH, H+)",
        "reaction2_reagents": "acetyl bromide (Base = LDA)",
    }
    options = {
        "A": {"A": "decahydro-7H-benzo[7]annulen-7-one", "B": "3-ethylpent-4-enoic acid"},
        "B": {"A": "decahydro-7H-benzo[7]annulen-7-one", "B": "lithium 3-ethylpent-4-enoate"},
        "C": {"A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one", "B": "lithium 3-ethylpent-4-enoate"},
        "D": {"A": "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one", "B": "3-ethylpent-4-enoic acid"}
    }
    smiles_map = {
        "1-vinylspiro[3.5]non-5-en-1-ol": "C=C[C]1(O)CC[C]12C=CCCC2",
        "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one": "C1CC2=C3C(CCC1)C(CC3)C(=O)C2",
        "decahydro-7H-benzo[7]annulen-7-one": "C1CCC2C(C1)CCCC(=O)CC2",
        "3-ethylpent-4-enoic acid": "C=CC(CC)CC(=O)O",
        "lithium 3-ethylpent-4-enoate": "C=CC(CC)CC(=O)[O-].[Li+]"
    }
    
    provided_answer_letter = "C"
    provided_answer_content = options[provided_answer_letter]
    errors = []

    # --- Constraint 1: Check Isomerization for Reaction 1 ---
    start_mat_1_name = "1-vinylspiro[3.5]non-5-en-1-ol"
    product_A_isomer_name = "(E)-bicyclo[5.3.1]undec-1(11)-en-4-one"
    product_A_non_isomer_name = "decahydro-7H-benzo[7]annulen-7-one"

    start_formula = get_molecular_formula(start_mat_1_name, smiles_map)
    isomer_formula = get_molecular_formula(product_A_isomer_name, smiles_map)
    non_isomer_formula = get_molecular_formula(product_A_non_isomer_name, smiles_map)

    if start_formula != isomer_formula:
        errors.append(f"Constraint 1 (Isomerization) is violated: The proposed product A '{product_A_isomer_name}' ({isomer_formula}) is not an isomer of the starting material '{start_mat_1_name}' ({start_formula}).")
    
    if start_formula == non_isomer_formula:
        errors.append(f"Constraint 1 (Isomerization) is violated: The reasoning to exclude '{product_A_non_isomer_name}' is flawed, as its formula ({non_isomer_formula}) matches the starting material ({start_formula}).")

    if provided_answer_content["A"] != product_A_isomer_name:
        errors.append(f"The final answer for product A is inconsistent with the isomerization constraint. It should be '{product_A_isomer_name}'.")

    # --- Constraint 2: Check Workup Conditions for Reaction 2 ---
    reaction_2_reagents = question["reaction2_reagents"]
    has_acid_workup = bool(re.search(r'\bH\+\b|H3O', reaction_2_reagents, re.IGNORECASE))

    product_B_name = provided_answer_content["B"]
    is_salt = "lithium" in product_B_name.lower() or "oate" in product_B_name.lower()

    if has_acid_workup:
        errors.append("Constraint 2 (Workup) is violated: The reasoning assumes no acidic workup, but the reagents for Reaction 2 seem to contain one.")
    elif not is_salt:
        errors.append(f"The final answer for product B is inconsistent with the workup conditions. Without an acidic workup, the product should be a salt, but the answer chose '{product_B_name}'.")

    # --- Final Verdict ---
    if not errors:
        # The logic holds and leads to the selected answer.
        return "Correct"
    else:
        # The logic is flawed or the answer doesn't follow the logic.
        error_message = "Incorrect. The following constraints or conditions were not satisfied:\n- "
        error_message += "\n- ".join(errors)
        return error_message

# Run the check and print the result
result = check_answer_correctness()
print(result)