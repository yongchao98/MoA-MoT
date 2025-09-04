def check_answer():
    """
    Checks the correctness of the answer by verifying chemical constraints.
    1. Checks if the molecular formula of the options matches C11H12O.
    2. Checks if the structure contains a p-tolyl group, which is necessary to form the product.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return "Skipping check: rdkit library is not installed. Please install it using 'pip install rdkit'."

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    # Define the options from the question. SMILES strings represent the structures.
    # The SMILES are generated from the IUPAC names.
    options = {
        "A": {"name": "2-styrylepoxide", 
              "smiles": "c1ccc(cc1)/C=C/C1OC1"},
        "B": {"name": "2-(4-methylstyryl)oxirane", 
              "smiles": "Cc1ccc(cc1)/C=C/C1OC1"},
        "C": {"name": "2-(1-phenylprop-1-en-2-yl)oxirane", 
              "smiles": "C=C(c1ccccc1)C1OC1"},
        "D": {"name": "2-methyl-3-styryloxirane", 
              "smiles": "Cc1oc1/C=C/c1ccccc1"}
    }

    # --- Define Constraints ---
    # Constraint 1: The molecular formula must be C11H12O.
    target_formula = "C11H12O"
    
    # Constraint 2: The structure must contain a p-tolyl group.
    # The product has a p-tolyl group, so the starting material must have one.
    # SMARTS pattern for a p-tolyl group (a methyl-substituted benzene ring with another attachment at the para position).
    p_tolyl_smarts = Chem.MolFromSmarts('Cc1ccc([!#1])cc1')

    # --- Verification Logic ---
    valid_candidates = []
    failure_reasons = {}

    for key, data in options.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            failure_reasons[key] = f"Option {key} ({data['name']}) has an invalid chemical structure representation (SMILES)."
            continue

        # Check Constraint 1: Molecular Formula
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        if formula != target_formula:
            failure_reasons[key] = f"Option {key} ({data['name']}) has the wrong molecular formula ({formula}). The required formula is {target_formula}."
            continue
            
        # Check Constraint 2: Structural Feature (p-tolyl group)
        if not mol.HasSubstructMatch(p_tolyl_smarts):
            failure_reasons[key] = f"Option {key} ({data['name']}) does not contain a p-tolyl group, which is required to form the product identified by the NMR spectra."
            continue

        # If a candidate passes all checks, add it to the list.
        valid_candidates.append(key)

    # --- Final Verdict ---
    if len(valid_candidates) == 1:
        correct_option = valid_candidates[0]
        if llm_answer == correct_option:
            return "Correct"
        else:
            reason = failure_reasons.get(llm_answer, f"Option {llm_answer} does not satisfy all constraints.")
            return f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option}. Reason: {reason}"
    elif len(valid_candidates) == 0:
        return "Incorrect. No option satisfies all the necessary chemical constraints."
    else:
        return f"Incorrect. The question is ambiguous as multiple options ({', '.join(valid_candidates)}) satisfy the constraints."

# Run the check and print the result.
result = check_answer()
print(result)