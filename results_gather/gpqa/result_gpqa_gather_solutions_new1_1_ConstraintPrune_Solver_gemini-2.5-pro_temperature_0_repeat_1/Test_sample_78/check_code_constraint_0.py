def check_correctness():
    """
    Checks the correctness of the answer by verifying the chemical logic.
    This function uses the RDKit library. To run it, you need to install RDKit:
    pip install rdkit
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return "Could not run the check because the 'rdkit' library is not installed. Please install it using 'pip install rdkit'."

    # --- Problem Constraints & Definitions ---
    # The starting material, Compound X, has the formula C11H12O.
    TARGET_FORMULA = "C11H12O"
    
    # The final answer given by the LLM.
    LLM_ANSWER = 'B'

    # Structures of the options for Compound X.
    options = {
        'A': {'name': '2-styrylepoxide', 'smiles': 'c1ccc(C=Cc2co2)cc1'},
        'B': {'name': '2-(4-methylstyryl)oxirane', 'smiles': 'Cc1ccc(C=Cc2co2)cc1'},
        'C': {'name': '2-methyl-3-styryloxirane', 'smiles': 'c1ccc(C=CC2C(C)O2)cc1'},
        'D': {'name': '2-(1-phenylprop-1-en-2-yl)oxirane', 'smiles': 'c1ccc(C=C(C)c2co2)cc1'}
    }

    # The product structure deduced from NMR data.
    product = {
        'name': '(E)-4-(p-tolyl)but-3-en-2-one',
        'smiles': 'CC(=O)/C=C/c1ccc(C)cc1'
    }

    # --- Step 1: Verify the product structure's key features ---
    # The product must be an isomer of the starting material and contain the p-tolyl group.
    product_mol = Chem.MolFromSmiles(product['smiles'])
    product_formula = Descriptors.rdMolDescriptors.CalcMolFormula(product_mol)
    
    if product_formula != TARGET_FORMULA:
        return f"Reason: The logic is flawed. The identified product {product['name']} has formula {product_formula}, which does not match the starting material's formula {TARGET_FORMULA}."

    # The p-tolyl group is the key structural feature identified from the NMR.
    p_tolyl_pattern = Chem.MolFromSmarts('c1(C)ccc(cc1)')
    if not product_mol.HasSubstructMatch(p_tolyl_pattern):
        return "Reason: The logic is flawed. The identified product does not contain a p-tolyl group, contradicting the NMR analysis."

    # --- Step 2: Filter the options based on the constraints ---
    valid_candidates = []
    for key, data in options.items():
        mol = Chem.MolFromSmiles(data['smiles'])
        if not mol:
            continue # Skip if SMILES is invalid

        # Constraint 1: Must have the correct molecular formula.
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
        formula_ok = (formula == TARGET_FORMULA)

        # Constraint 2: Must contain the p-tolyl group to form the product.
        has_p_tolyl = mol.HasSubstructMatch(p_tolyl_pattern)

        if formula_ok and has_p_tolyl:
            valid_candidates.append(key)

    # --- Step 3: Draw a conclusion ---
    if len(valid_candidates) != 1:
        return f"Reason: The filtering logic is inconclusive. Found {len(valid_candidates)} candidates ({valid_candidates}) that match the criteria, but expected exactly one."

    derived_correct_answer = valid_candidates[0]

    if LLM_ANSWER == derived_correct_answer:
        return "Correct"
    else:
        # Explain why the LLM's answer is wrong by checking which constraint it fails.
        failed_mol = Chem.MolFromSmiles(options[LLM_ANSWER]['smiles'])
        failed_formula = Descriptors.rdMolDescriptors.CalcMolFormula(failed_mol)
        if failed_formula != TARGET_FORMULA:
            return f"Reason: The answer '{LLM_ANSWER}' is incorrect. Its molecular formula is {failed_formula}, but the question specifies {TARGET_FORMULA}."
        
        if not failed_mol.HasSubstructMatch(p_tolyl_pattern):
            return f"Reason: The answer '{LLM_ANSWER}' is incorrect. It does not contain a p-tolyl group, which is required to form the product identified from the NMR data."
            
        return f"Reason: The answer '{LLM_ANSWER}' is incorrect. The correct answer is '{derived_correct_answer}'."

# Execute the check and print the result.
print(check_correctness())