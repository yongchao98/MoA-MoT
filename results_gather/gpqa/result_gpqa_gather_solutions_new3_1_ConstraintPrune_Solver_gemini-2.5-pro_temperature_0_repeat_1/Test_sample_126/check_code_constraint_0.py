# Install rdkit if you don't have it:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement_product():
    """
    Checks the correctness of the proposed answer for the Cope rearrangement of 5-butylnona-2,6-diene.
    """
    # Define the reactant and options using their IUPAC names and SMILES strings
    # SMILES strings are a standard way to represent chemical structures.
    reactant = {
        "name": "5-butylnona-2,6-diene",
        "smiles": "CCC=CC(CCCC)C=CCC"
    }
    
    options = {
        "A": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"},
        "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"} # Same as B
    }
    
    proposed_answer_key = "A"
    
    # --- Step 1: Isomer Check ---
    print("--- Constraint 1: Product must be an isomer of the reactant. ---")
    
    try:
        reactant_mol = Chem.MolFromSmiles(reactant["smiles"])
        if not reactant_mol:
            return f"Error: Could not parse reactant SMILES: {reactant['smiles']}"
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)
        print(f"Reactant: {reactant['name']} has formula {reactant_formula}")

        correct_isomer = True
        for key, data in options.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if not mol:
                return f"Error: Could not parse SMILES for option {key}: {data['smiles']}"
            
            formula = rdMolDescriptors.CalcMolFormula(mol)
            is_isomer = (formula == reactant_formula)
            
            print(f"Option {key}: {data['name']} has formula {formula}. Isomer? -> {is_isomer}")
            
            if key == proposed_answer_key and not is_isomer:
                return f"Incorrect: The proposed answer {proposed_answer_key} is not an isomer of the reactant."
            
            if not is_isomer and key != 'C': # We expect C to be a non-isomer
                 correct_isomer = False

        if not correct_isomer:
             print("Warning: An unexpected non-isomer was found.")
        print("Isomer Check Result: PASSED. Proposed answer A is a valid isomer. Option C is correctly identified as a non-isomer.\n")

    except Exception as e:
        return f"An error occurred during the isomer check: {e}"

    # --- Step 2: Mechanism Plausibility Check ---
    print("--- Constraint 2: Product must be a 1,5-diene, as predicted by the Cope rearrangement mechanism. ---")
    
    # A 1,5-diene has the pattern C=C-C-C-C=C
    # A 1,6-diene has the pattern C=C-C-C-C-C=C
    # We use SMARTS patterns to find these substructures.
    pattern_1_5_diene = Chem.MolFromSmarts('[#6]=[#6]-[#6]-[#6]-[#6]=[#6]')
    pattern_1_6_diene = Chem.MolFromSmarts('[#6]=[#6]-[#6]-[#6]-[#6]-[#6]=[#6]')

    answer_satisfies_mechanism = False
    
    for key, data in options.items():
        # Skip the non-isomer
        if key == 'C':
            continue
            
        mol = Chem.MolFromSmiles(data["smiles"])
        is_1_5_diene = mol.HasSubstructMatch(pattern_1_5_diene)
        is_1_6_diene = mol.HasSubstructMatch(pattern_1_6_diene)
        
        diene_type = "1,5-diene" if is_1_5_diene else "1,6-diene" if is_1_6_diene else "Other"
        print(f"Option {key} ({data['name']}) is a: {diene_type}")
        
        if key == proposed_answer_key:
            if is_1_5_diene:
                answer_satisfies_mechanism = True
            else:
                return f"Incorrect: The proposed answer {proposed_answer_key} is a {diene_type}, but the Cope rearrangement should produce a 1,5-diene."

    if answer_satisfies_mechanism:
        print("Mechanism Plausibility Check Result: PASSED. Proposed answer A is a 1,5-diene, consistent with the mechanism.\n")
    else:
        # This case should not be reached if the logic is correct
        return "Incorrect: The proposed answer did not satisfy the mechanism check."

    return "Correct"

# Run the check
result = check_cope_rearrangement_product()
print(f"\nFinal Verdict: {result}")

# Final output based on the script's logic
if result == "Correct":
    print("The provided answer <<<A>>> is correct because it is an isomer of the reactant and its structure (a 1,5-diene) is consistent with the product of a Cope rearrangement.")
else:
    print(f"The provided answer is incorrect. Reason: {result}")
