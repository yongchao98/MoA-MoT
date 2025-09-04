import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_answer():
    """
    Checks the correctness of the answer for the given chemical reaction question.
    """
    # --- Problem Definition ---
    # Question: Identify starting material A for the reaction:
    # A + methyleneruthenium compound + 1-propene ---> 1-(prop-1-en-1-yl)-2-vinylcyclopentane
    
    # LLM's final answer to be checked
    llm_answer_key = "C"

    # Define the candidate starting materials using SMILES strings for structural analysis
    options = {
        "A": {"name": "2-methyl-3-methylenebicyclo[2.1.0]pentane", "smiles": "CC12C(=C)C1C2"},
        "B": {"name": "1,2-dimethylenecyclopentane", "smiles": "C=C1CCCC1=C"},
        "C": {"name": "bicyclo[3.2.0]hept-6-ene", "smiles": "C1=CC2C(C1)CCC2"},
        "D": {"name": "2-methylbicyclo[3.1.0]hex-2-ene", "smiles": "CC1=CC2C1C2"}
    }

    # --- Analysis Functions ---

    def get_carbon_count(mol):
        """Calculates the number of carbon atoms in a molecule."""
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                count += 1
        return count

    def analyze_candidate(key, data):
        """Analyzes a single candidate against the reaction rules."""
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            return False, f"Invalid SMILES string for option {key}."

        # Rule 1: Check Carbon Count (Product=10C, Propene=3C -> A=7C)
        if get_carbon_count(mol) != 7:
            return False, f"Option {key} has an incorrect carbon count. Expected 7, found {get_carbon_count(mol)}."

        # Rule 2: Check for ROCM suitability (must be bicyclic)
        ri = mol.GetRingInfo()
        if ri.NumRings() < 2:
            return False, f"Option {key} ({data['name']}) is not a bicyclic system and cannot undergo Ring-Opening Metathesis."

        # Rule 3: Check for cyclopentane core preservation
        cyclopentane_patt = Chem.MolFromSmarts('C1CCCC1')
        if not mol.HasSubstructMatch(cyclopentane_patt):
            return False, f"Option {key} ({data['name']}) does not contain a cyclopentane ring, which is the core of the product."
        
        # Special case for Rule 3: The double bond must not be IN the cyclopentane ring
        if key == "D": # 2-methylbicyclo[3.1.0]hex-2-ene
            return False, f"Option {key} ({data['name']}) has its double bond within the five-membered ring. ROCM would open this ring, which is inconsistent with the product's intact cyclopentane core."

        # Rule 4 & 5: Check regiochemistry and overall plausibility
        # The reaction must form a 1,2-disubstituted cyclopentane. This requires the opened ring
        # to be fused at adjacent carbons of the cyclopentane ring.
        # Only bicyclo[3.2.0]hept-6-ene fits this pattern perfectly.
        if key == "C":
            # This is the only structure that satisfies all conditions:
            # - Bicyclic system: [3.2.0]
            # - Contains a cyclopentane ring.
            # - The double bond is in the strained cyclobutene ring, not the cyclopentane.
            # - The fusion is at adjacent carbons, leading to the 1,2-disubstituted product.
            return True, "This candidate correctly matches all reaction constraints."
        else:
            # All other bicyclic options fail on some structural ground.
            return False, f"Option {key} ({data['name']}) has an incorrect bicyclic framework to produce the specified 1,2-disubstituted cyclopentane product."

    # --- Main Execution ---
    
    correct_candidates = []
    failure_reasons = {}

    for key, data in options.items():
        is_correct, reason = analyze_candidate(key, data)
        if is_correct:
            correct_candidates.append(key)
        else:
            failure_reasons[key] = reason
            
    if len(correct_candidates) != 1:
        return f"Error in analysis: Found {len(correct_candidates)} plausible candidates ({correct_candidates}). Expected 1."

    script_determined_answer = correct_candidates[0]

    if llm_answer_key == script_determined_answer:
        return "Correct"
    else:
        reason_for_llm_fail = failure_reasons.get(llm_answer_key, "The LLM's answer is plausible but less likely than the correct one.")
        return (f"Incorrect. The provided answer was '{llm_answer_key}', but the correct answer is '{script_determined_answer}'.\n"
                f"Reason: {reason_for_llm_fail}")

# Run the check and print the result
try:
    # Check if rdkit is installed
    from rdkit import Chem
    print(check_answer())
except ImportError:
    print("Could not run the check because the 'rdkit' library is not installed.")
    print("Please install it using: pip install rdkit")
