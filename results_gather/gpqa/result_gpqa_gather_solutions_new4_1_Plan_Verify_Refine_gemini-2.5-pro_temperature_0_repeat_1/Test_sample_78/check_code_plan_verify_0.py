from rdkit import Chem
from rdkit.Chem import Descriptors

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by:
    1. Verifying the logical deduction about the product's structure.
    2. Applying the key logical constraints (molecular formula, conserved substructure) to the provided options.
    3. Comparing the logically derived correct option with the LLM's answer.
    """
    
    # --- Problem Constraints & LLM Answer ---
    target_formula = "C11H12O"
    llm_provided_answer = "B"

    # --- Structure Definitions (SMILES and SMARTS) ---
    # Product deduced from NMR: (E)-4-(4-methylphenyl)but-3-en-2-one
    product_smiles = "CC(=O)/C=C/c1ccc(C)cc1"
    
    # Options for Compound X (interpreted from IUPAC names)
    options_smiles = {
        "A": "c1ccc(cc1)/C=C/C2OC2C",          # 2-methyl-3-styryloxirane
        "B": "Cc1ccc(cc1)/C=C/C2CO2",        # 2-(4-methylstyryl)oxirane
        "C": "c1ccc(cc1)/C=C/C2CO2",        # 2-styrylepoxide (a.k.a. 2-styryloxirane)
        "D": "CC(=Cc1ccccc1)C1OC1"           # A plausible C12H12O structure for 2-(1-phenylprop-1-en-2-yl)oxirane
    }
    
    # Substructure to be checked (p-tolyl group)
    p_tolyl_smarts = "[#6]c1ccc(C)cc1"

    # --- Step 1: Verify the deduced product structure is consistent ---
    product_mol = Chem.MolFromSmiles(product_smiles)
    if Chem.rdMolDescriptors.CalcMolFormula(product_mol) != target_formula:
        return "Internal Check Error: The reference product structure used for checking has an incorrect formula."
    if not product_mol.HasSubstructMatch(Chem.MolFromSmarts(p_tolyl_smarts)):
        return "Internal Check Error: The reference product structure is missing the p-tolyl group, which contradicts the NMR analysis."
    
    # --- Step 2: Identify the correct starting material based on logical constraints ---
    correct_candidate = None
    reasons_for_rejection = {}

    for key, smiles in options_smiles.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            reasons_for_rejection[key] = f"Could not parse the chemical structure (SMILES: {smiles})."
            continue

        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        has_p_tolyl = mol.HasSubstructMatch(Chem.MolFromSmarts(p_tolyl_smarts))

        # Constraint 1: Must be an isomer of the product (have the target formula)
        if formula != target_formula:
            reasons_for_rejection[key] = f"Has incorrect molecular formula ({formula}). The question specifies the formula is {target_formula}."
            continue
        
        # Constraint 2: Must contain the p-tolyl group (conserved from product)
        if not has_p_tolyl:
            reasons_for_rejection[key] = "Does not contain the required p-tolyl group, which is clearly present in the product's NMR spectrum."
            continue
        
        # If all constraints are met, this is our candidate
        if correct_candidate is not None:
            # This case would indicate an ambiguous question
            return "Error: The question is ambiguous as multiple options fit the logical criteria."
        correct_candidate = key

    # --- Step 3: Compare the identified correct candidate with the LLM's answer ---
    if correct_candidate is None:
        return "Error: None of the options fit the logical criteria to be the correct starting material."

    if llm_provided_answer == correct_candidate:
        return "Correct"
    else:
        rejection_reason = reasons_for_rejection.get(llm_provided_answer, "It does not meet the required structural constraints.")
        return (f"Incorrect. The provided answer is {llm_provided_answer}, but the only logical choice is {correct_candidate}. "
                f"Reason: Option {correct_candidate} is the only candidate that has the correct molecular formula ({target_formula}) and contains the p-tolyl group found in the product. "
                f"The chosen option {llm_provided_answer} is wrong because: {rejection_reason}")

# Run the check
result = check_correctness()
print(result)