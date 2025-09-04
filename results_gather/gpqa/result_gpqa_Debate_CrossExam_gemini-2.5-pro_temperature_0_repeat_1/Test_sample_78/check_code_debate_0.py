import sys
from io import StringIO

# Use a try-except block to handle the case where rdkit is not installed.
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'")
    # If rdkit is not available, we cannot perform the check.
    # We will return a message indicating the dependency is missing.
    # In a real scenario, this would cause the check to fail.
    # For this demonstration, we'll assume it's a setup issue.
    sys.exit("Missing dependency: rdkit")

def calculate_molecular_formula(mol):
    """Calculates the molecular formula from an RDKit molecule object."""
    return rdMolDescriptors.CalcMolFormula(mol)

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying its claims
    about molecular structures and formulas.
    """
    # --- Step 1: Define the problem's data and the answer's claims ---

    # Given data from the question
    target_formula_str = "C11H12O"
    
    # Options for Compound X (starting material) represented by SMILES strings
    # SMILES strings are a standard way to represent chemical structures.
    options = {
        'A': {'name': '2-methyl-3-styryloxirane', 'smiles': 'CC1OC1/C=C/c2ccccc2'},
        'B': {'name': '2-(4-methylstyryl)oxirane', 'smiles': 'C1OC1/C=C/c2ccc(C)cc2'},
        'C': {'name': '2-(1-phenylprop-1-en-2-yl)oxirane', 'smiles': 'C1OC1C(=C/c2ccccc2)C'},
        'D': {'name': '2-styrylepoxide', 'smiles': 'c1ccc(/C=C/C2CO2)cc1'}
    }

    # The product structure identified in the answer
    product_smiles = 'CC(=O)/C=C/c1ccc(C)cc1'
    product_name = '(E)-4-(p-tolyl)but-3-en-2-one'
    
    # The final answer provided by the LLM
    llm_answer = 'B'

    # --- Step 2: Verify the analysis of the options' molecular formulas ---
    
    # The answer claims A, B, C have formula C11H12O and D has C10H10O.
    # Let's check this.
    for key, data in options.items():
        mol = Chem.MolFromSmiles(data['smiles'])
        formula = calculate_molecular_formula(mol)
        
        if key in ['A', 'B', 'C']:
            if formula != target_formula_str:
                return f"Incorrect: The answer assumes option {key} has the formula {target_formula_str}, but it is actually {formula}."
        elif key == 'D':
            if formula == target_formula_str:
                return f"Incorrect: The answer correctly dismisses option {key} for having the wrong formula, but the calculated formula {formula} matches the target."
            # The answer states D is C10H10O, which our SMILES confirms.
            if formula != "C10H10O":
                 return f"Incorrect: The answer claims option D has formula C10H10O, but the code calculates {formula}."

    # --- Step 3: Verify the analysis of the product structure ---

    # The answer claims the product is (E)-4-(p-tolyl)but-3-en-2-one and that it matches the NMR.
    # Let's verify its structural features.
    product_mol = Chem.MolFromSmiles(product_smiles)
    
    # 3a. Check product's molecular formula
    product_formula = calculate_molecular_formula(product_mol)
    if product_formula != target_formula_str:
        return f"Incorrect: The proposed product '{product_name}' has formula {product_formula}, which does not match the starting material's formula {target_formula_str}."

    # 3b. Check for key structural features that explain the NMR data.
    # p-tolyl group: a benzene ring with a methyl group. SMARTS: c1ccc(C)cc1
    p_tolyl_pattern = Chem.MolFromSmarts('c1ccc(C)cc1')
    if not product_mol.HasSubstructMatch(p_tolyl_pattern):
        return f"Incorrect: The proposed product '{product_name}' does not contain a p-tolyl group, which contradicts the NMR interpretation of a para-substituted ring and a methyl singlet."

    # Acetyl group: CH3-C=O. SMARTS: [CH3]C=O
    acetyl_pattern = Chem.MolFromSmarts('[CH3]C=O')
    if not product_mol.HasSubstructMatch(acetyl_pattern):
        return f"Incorrect: The proposed product '{product_name}' does not contain an acetyl group, which contradicts the NMR interpretation of a ketone (13C) and a methyl singlet (1H)."

    # The presence of two distinct methyl groups (one tolyl, one acetyl) matches the two 1H NMR singlets.
    # The presence of the p-tolyl group matches the two 2H doublets.
    # The presence of the trans C=C bond between the ring and carbonyl matches the two 1H doublets.
    # The structural analysis in the answer is consistent with the data.

    # --- Step 4: Verify the logical link between reactant and product ---

    # The core argument: The reaction is an isomerization, so the p-tolyl group in the product
    # must have been present in the starting material.
    # Let's check which of the valid options (A, B, C) contains this group.
    
    correct_reactant_found = False
    for key, data in options.items():
        # Skip option D as it has the wrong formula
        if key == 'D':
            continue
            
        mol = Chem.MolFromSmiles(data['smiles'])
        has_p_tolyl = mol.HasSubstructMatch(p_tolyl_pattern)
        
        if key == llm_answer:
            if not has_p_tolyl:
                return f"Incorrect: The answer selects option {key}, but it does not contain the p-tolyl group required to form the product."
            else:
                correct_reactant_found = True
        else: # For other options
            if has_p_tolyl:
                return f"Incorrect: The answer dismisses option {key}, but it also contains a p-tolyl group, making it a potential candidate based on the answer's own logic."

    if not correct_reactant_found:
        return f"Incorrect: The logic failed to identify the chosen answer {llm_answer} as the correct precursor."

    # --- Step 5: Final Conclusion ---
    # All checks passed. The answer's reasoning is sound and verifiable.
    # 1. It correctly analyzes the molecular formulas of the options.
    # 2. It correctly deduces the product structure from the NMR data.
    # 3. It correctly applies chemical logic (conservation of the p-tolyl group) to link the product back to the starting material.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)
