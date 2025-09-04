def check_chemistry_answer():
    """
    Checks the correctness of the answer to the chemistry question by verifying
    molecular formulas and substructures of the given options.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return "Cannot perform check: The 'rdkit' library is not installed. Please install it using 'pip install rdkit'."

    # --- Problem Definition ---
    
    # The final answer provided for checking
    final_answer_choice = 'B'
    
    # Constraints derived from the question
    required_formula = "C11H12O"
    # The p-tolyl substructure is identified from the product's NMR data.
    # We use a SMARTS string to represent this substructure.
    p_tolyl_smarts = "c1(C)ccc(cc1)"
    p_tolyl_pattern = Chem.MolFromSmarts(p_tolyl_smarts)

    # Define the options with their chemical names and SMILES representations
    options = {
        'A': {
            "name": "2-styrylepoxide",
            "smiles": "c1ccc(cc1)C=CC2CO2"  # 2-(phenylethenyl)oxirane
        },
        'B': {
            "name": "2-(4-methylstyryl)oxirane",
            "smiles": "Cc1ccc(cc1)C=CC2CO2"  # 2-[2-(4-methylphenyl)ethenyl]oxirane
        },
        'C': {
            "name": "2-(1-phenylprop-1-en-2-yl)oxirane",
            "smiles": "c1ccccc1C=C(C)C2CO2"
        },
        'D': {
            "name": "2-methyl-3-styryloxirane",
            "smiles": "CC1C(C=Cc2ccccc2)O1"
        }
    }

    # --- Verification Logic ---
    
    correctly_identified_options = []
    analysis_log = []

    for key, data in options.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            analysis_log.append(f"Option {key} ('{data['name']}') has an invalid SMILES string.")
            continue

        # Constraint 1: Check the molecular formula
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        has_correct_formula = (formula == required_formula)

        # Constraint 2: Check for the p-tolyl substructure
        has_p_tolyl_group = mol.HasSubstructMatch(p_tolyl_pattern)
        
        analysis_log.append(
            f"  - Option {key} ({data['name']}): Formula is {formula} (Correct: {has_correct_formula}). "
            f"Contains p-tolyl group: {has_p_tolyl_group}."
        )

        # Check if this option satisfies all conditions
        if has_correct_formula and has_p_tolyl_group:
            correctly_identified_options.append(key)

    # --- Final Verdict ---
    
    # There should be exactly one option that satisfies all constraints
    if len(correctly_identified_options) == 1:
        derived_correct_choice = correctly_identified_options[0]
        if derived_correct_choice == final_answer_choice:
            return "Correct"
        else:
            reason = (f"The provided answer '{final_answer_choice}' is incorrect. "
                      f"The analysis shows that the correct answer is '{derived_correct_choice}'.\n"
                      f"Reasoning:\n"
                      f"1. The starting material must have the molecular formula {required_formula}.\n"
                      f"2. The product's NMR data indicates a p-tolyl group, which must also be present in the starting material.\n"
                      f"Analysis Details:\n" + "\n".join(analysis_log))
            return reason
    elif len(correctly_identified_options) == 0:
        reason = ("No option satisfies all the required constraints.\n"
                  "Analysis Details:\n" + "\n".join(analysis_log))
        return f"Incorrect. {reason}"
    else: # More than one option satisfies the constraints
        reason = (f"Multiple options ({', '.join(correctly_identified_options)}) satisfy the constraints, "
                  "indicating an ambiguity in the problem or analysis.\n"
                  "Analysis Details:\n" + "\n".join(analysis_log))
        return f"Incorrect. {reason}"

# Run the check and print the result
print(check_chemistry_answer())