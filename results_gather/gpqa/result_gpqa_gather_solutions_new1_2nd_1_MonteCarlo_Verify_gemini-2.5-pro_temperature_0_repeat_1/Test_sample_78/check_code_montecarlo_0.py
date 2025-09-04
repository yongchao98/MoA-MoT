# You may need to install rdkit: pip install rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def check_correctness():
    """
    Checks the correctness of the provided answer by verifying the chemical constraints.
    """
    # The final answer from the LLM to be checked
    llm_answer = "D"

    # --- Step 1: Define the problem constraints and options ---

    # Constraint 1: Molecular formula of the starting material (Compound X)
    target_formula = "C11H12O"

    # Constraint 2: The product's NMR shows a p-tolyl group. This group must be
    # conserved from the starting material. We define a SMARTS pattern to find it.
    # This pattern looks for a benzene ring with a methyl group, where the ring is
    # attached to another non-hydrogen atom at the para position.
    ptolyl_smarts = "[#6]c1ccc(C)cc1"
    ptolyl_mol_pattern = Chem.MolFromSmarts(ptolyl_smarts)

    # Define the options from the question with their chemical structures as SMILES strings
    options = {
        "A": {
            "name": "2-(1-phenylprop-1-en-2-yl)oxirane",
            "smiles": "CC(=Cc1ccccc1)C2OC2"
        },
        "B": {
            "name": "2-methyl-3-styryloxirane",
            "smiles": "CC1OC1C=Cc1ccccc1"
        },
        "C": {
            "name": "2-styrylepoxide",
            "smiles": "c1ccc(C=CC2OC2)cc1"
        },
        "D": {
            "name": "2-(4-methylstyryl)oxirane",
            "smiles": "Cc1ccc(C=CC2OC2)cc1"
        }
    }

    # --- Step 2: Verify each option against the constraints ---
    
    valid_candidates = []
    reasons_for_failure = {}

    for letter, data in options.items():
        try:
            mol = Chem.MolFromSmiles(data["smiles"])
            if mol is None:
                reasons_for_failure[letter] = f"Option {letter} ('{data['name']}') has an invalid chemical structure representation (SMILES)."
                continue

            # Check 1: Molecular Formula
            formula = CalcMolFormula(mol)
            if formula != target_formula:
                reasons_for_failure[letter] = f"Option {letter} ('{data['name']}') has the wrong molecular formula ({formula}). The required formula is {target_formula}."
                continue

            # Check 2: Presence of the required p-tolyl group
            if not mol.HasSubstructMatch(ptolyl_mol_pattern):
                reasons_for_failure[letter] = f"Option {letter} ('{data['name']}') has the correct formula but lacks the required p-tolyl group. It has a phenyl group instead."
                continue
            
            # If a candidate passes all checks, it is valid.
            valid_candidates.append(letter)

        except Exception as e:
            reasons_for_failure[letter] = f"An error occurred while processing option {letter}: {e}"

    # --- Step 3: Determine the correct answer and compare with the LLM's answer ---

    if len(valid_candidates) == 1:
        correct_answer = valid_candidates[0]
        if llm_answer == correct_answer:
            return "Correct"
        else:
            failure_reason = reasons_for_failure.get(llm_answer, f"The provided answer {llm_answer} does not meet the problem's constraints.")
            return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_answer}.\n"
                    f"Reason: {failure_reason}\n"
                    f"Only option {correct_answer} ('{options[correct_answer]['name']}') satisfies all conditions: it has the molecular formula {target_formula} and contains the necessary p-tolyl group to form the product identified by the NMR spectra.")
    elif len(valid_candidates) == 0:
        return f"Incorrect. No option satisfies all the required conditions. The provided answer was {llm_answer}. Analysis of options: {reasons_for_failure}"
    else: # More than one valid candidate
        return f"Incorrect. The problem is ambiguous as multiple options ({', '.join(valid_candidates)}) satisfy the conditions. The provided answer was {llm_answer}."

# Run the check
result = check_correctness()
print(result)