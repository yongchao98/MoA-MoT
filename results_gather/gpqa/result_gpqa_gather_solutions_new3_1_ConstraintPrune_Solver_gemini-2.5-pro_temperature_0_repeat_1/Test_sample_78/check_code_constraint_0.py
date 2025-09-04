def check_chemistry_answer():
    """
    This function checks the correctness of the LLM's answer by:
    1. Defining the constraints from the problem (molecular formula, required substructure).
    2. Representing the candidate molecules using SMILES strings.
    3. Using the RDKit library to verify which candidate satisfies all constraints.
    4. Comparing the valid candidate with the LLM's provided answer.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit-pypi'"

    # --- Step 1: Define Constraints from the Problem ---

    # The molecular formula of Compound X is given as C11H12O.
    TARGET_FORMULA = "C11H12O"

    # The NMR analysis of the product clearly indicates a p-tolyl group.
    # Since the reaction is an isomerization, the starting material must also have this group.
    # We define this required substructure using a SMARTS string.
    # SMARTS for a p-tolyl group: a methyl-substituted benzene ring, para to an attachment point (*).
    P_TOLYL_SMARTS = "c1(C)ccc(c*)cc1"

    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'D'

    # --- Step 2: Define the Candidate Molecules ---

    # The options are represented by their SMILES strings.
    options = {
        'A': "CC(=Cc1ccccc1)C1CO1",  # 2-(1-phenylprop-1-en-2-yl)oxirane
        'B': "c1ccc(C=CC2CO2)cc1",    # 2-styrylepoxide
        'C': "CC1OC1C=Cc1ccccc1",    # 2-methyl-3-styryloxirane
        'D': "Cc1ccc(C=CC2CO2)cc1"     # 2-(4-methylstyryl)oxirane
    }

    # --- Step 3: Validate Each Option Against Constraints ---

    valid_options = []
    rejection_reasons = {}

    for key, smiles in options.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            rejection_reasons[key] = "Invalid chemical structure (SMILES)."
            continue

        # Constraint 1: Check the molecular formula.
        formula = rdMolDescriptors.CalcMolFormula(mol)
        if formula != TARGET_FORMULA:
            rejection_reasons[key] = f"Incorrect molecular formula. Expected {TARGET_FORMULA}, but got {formula}."
            continue

        # Constraint 2: Check for the presence of the p-tolyl group.
        p_tolyl_pattern = Chem.MolFromSmarts(P_TOLYL_SMARTS)
        if not mol.HasSubstructMatch(p_tolyl_pattern):
            rejection_reasons[key] = "Missing required substructure: Does not contain a p-tolyl group."
            continue

        # If both constraints are passed, the option is valid.
        valid_options.append(key)

    # --- Step 4: Check the LLM's Answer ---

    # A single, unique option should be valid.
    if len(valid_options) != 1:
        return f"Error in problem analysis: Found {len(valid_options)} options ({valid_options}) that satisfy the constraints. The problem should have a unique answer."

    correct_answer = valid_options[0]

    if llm_final_answer == correct_answer:
        return "Correct"
    else:
        reason = rejection_reasons.get(llm_final_answer, "It does not satisfy the problem's constraints.")
        return f"Incorrect. The provided answer '{llm_final_answer}' is wrong. Reason: {reason}. The correct answer is '{correct_answer}' because it is the only option with the correct molecular formula and the required p-tolyl group."

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)