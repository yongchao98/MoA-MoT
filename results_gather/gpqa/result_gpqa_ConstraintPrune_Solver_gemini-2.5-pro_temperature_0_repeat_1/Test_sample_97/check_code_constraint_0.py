import sys

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer for the reaction:
    A + 1-propene -> 1-(prop-1-en-1-yl)-2-vinylcyclopentane

    The function verifies the step-by-step pruning logic presented in the answer
    using the RDKit chemistry library.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return ("CRITICAL_ERROR: The 'rdkit' library is required to perform this check. "
                "Please install it (e.g., 'pip install rdkit-pypi') and run the script again.")

    # --- Data Definition ---
    # Define molecules using SMILES for unambiguous chemical representation.
    candidates = {
        "A": {"name": "2-methylbicyclo[3.1.0]hex-2-ene", "smiles": "CC1=CC2C1C2"},
        "B": {"name": "bicyclo[3.2.0]hept-6-ene", "smiles": "C1CC2C(C1)C=C2"},
        "C": {"name": "1,2-dimethylenecyclopentane", "smiles": "C=C1CCCC1=C"},
        "D": {"name": "2-methyl-3-methylenebicyclo[2.1.0]pentane", "smiles": "C=C1C(C)C2C1C2"}
    }
    product_smiles = "C=CC1CCCC1C=CC"
    propene_smiles = "C=CC"
    correct_answer_key = "B"

    # --- Helper Functions ---
    def get_carbon_count(mol):
        return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    def is_bicyclic(mol):
        return mol.GetRingInfo().NumRings() == 2

    def has_cyclopentane_scaffold(mol):
        return any(len(ring) == 5 for ring in mol.GetRingInfo().AtomRings())

    def is_double_bond_in_cyclopentane(mol):
        double_bond_match = mol.GetSubstructMatch(Chem.MolFromSmarts("C=C"))
        if not double_bond_match: return False
        for ring in mol.GetRingInfo().AtomRings():
            if len(ring) == 5 and double_bond_match[0] in ring and double_bond_match[1] in ring:
                return True
        return False

    # --- Verification Logic ---
    # 1. Parse all molecules and check for errors
    try:
        product_mol = Chem.MolFromSmiles(product_smiles)
        propene_mol = Chem.MolFromSmiles(propene_smiles)
        for key in candidates:
            candidates[key]['mol'] = Chem.MolFromSmiles(candidates[key]['smiles'])
            if candidates[key]['mol'] is None:
                return f"Internal Error: Failed to parse SMILES for candidate {key}."
    except Exception:
        return "Internal Error: An exception occurred during molecule parsing."

    # Constraint 1: Carbon Count. A must have 7 carbons.
    required_carbons = get_carbon_count(product_mol) - get_carbon_count(propene_mol)
    if required_carbons != 7:
        return "Incorrect analysis: The carbon count for starting material 'A' should be 7, but calculation yields a different number."
    for key, data in candidates.items():
        if get_carbon_count(data['mol']) != 7:
            return f"Incorrect analysis: Candidate {key} does not have 7 carbons as stated in the provided answer."

    # Constraint 2: Reaction Type (ROCM). Prune non-bicyclic candidates.
    passing_candidates_c2 = {k: v for k, v in candidates.items() if is_bicyclic(v['mol'])}
    if 'C' in passing_candidates_c2 or 'A' not in passing_candidates_c2 or 'B' not in passing_candidates_c2 or 'D' not in passing_candidates_c2:
        return "Incorrect analysis for Constraint 2: The check for bicyclic structures is flawed. Only C should be pruned at this step."

    # Constraint 3: Product Core Structure. Prune candidates that don't yield a cyclopentane core.
    # This means pruning those without a cyclopentane ring or where the double bond is in the cyclopentane ring.
    passing_candidates_c3 = {}
    for key, data in passing_candidates_c2.items():
        if has_cyclopentane_scaffold(data['mol']) and not is_double_bond_in_cyclopentane(data['mol']):
            passing_candidates_c3[key] = data
    
    if 'A' in passing_candidates_c3 or 'D' in passing_candidates_c3 or 'B' not in passing_candidates_c3:
        return "Incorrect analysis for Constraint 3: The logic for preserving the cyclopentane core is flawed. Candidate A (DB in ring) and D (no 5-ring) should be pruned, leaving only B."

    # Final check on the remaining candidate
    if len(passing_candidates_c3) != 1 or correct_answer_key not in passing_candidates_c3:
        return f"Incorrect conclusion: After applying all constraints, the remaining candidate is not {correct_answer_key} as stated. Remaining candidates: {list(passing_candidates_c3.keys())}"

    # If all checks pass, the reasoning is sound and the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_answer_correctness()
print(result)