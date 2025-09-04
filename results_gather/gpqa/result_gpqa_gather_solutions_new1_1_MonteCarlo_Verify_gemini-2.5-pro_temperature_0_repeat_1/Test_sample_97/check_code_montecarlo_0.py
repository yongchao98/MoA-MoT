try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit'")
    # Create a dummy class to avoid crashing the script if rdkit is not present.
    class Chem:
        @staticmethod
        def MolFromSmiles(smiles): return None
        @staticmethod
        def GetSSSR(mol): return 0
        class rdchem:
            class BondType:
                DOUBLE = None

def check_correctness():
    """
    Checks the correctness of the LLM's answer by applying chemical constraints.
    """

    # --- Data Setup ---
    # The question identifies the starting material 'A' for the reaction:
    # A + 1-propene --(Ru catalyst)--> 1-(prop-1-en-1-yl)-2-vinylcyclopentane
    options = {
        "A": {"name": "2-methylbicyclo[3.1.0]hex-2-ene", "smiles": "CC1=CC2C1C2"},
        "B": {"name": "2-methyl-3-methylenebicyclo[2.1.0]pentane", "smiles": "C=C1C(C)C2C1C2"},
        "C": {"name": "1,2-dimethylenecyclopentane", "smiles": "C=C1CCCC1=C"},
        "D": {"name": "bicyclo[3.2.0]hept-6-ene", "smiles": "C1=CC2C(C1)CC2"}
    }
    llm_answer_key = "D"

    # --- Constraint Checking Function ---
    def check_candidate(smiles: str) -> (bool, str):
        """
        Checks a single candidate molecule against the reaction constraints.
        Returns a tuple: (passes_checks, reason_for_failure_or_success).
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False, "Invalid chemical structure (SMILES)."

        # Constraint 1: Carbon Count must be 7.
        # Product (C10H16) - 1-propene (C3H6) = A (C7H10)
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if num_carbons != 7:
            return False, f"Fails carbon count. Expected 7, found {num_carbons}."

        # Constraint 2: Must be a bicyclic alkene for Ring-Opening Cross-Metathesis (ROCM).
        if Chem.GetSSSR(mol) != 2:
            return False, "Fails structural check: Not a bicyclic compound."

        # Constraint 3: Must contain a 5-membered ring to form the cyclopentane core.
        ring_info = mol.GetRingInfo()
        ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
        if 5 not in ring_sizes:
            return False, "Fails structural check: Does not contain a 5-membered ring."

        # Constraint 4: The double bond must be outside the 5-membered ring for it to survive.
        double_bonds = [b for b in mol.GetBonds() if b.GetBondType() == Chem.rdchem.BondType.DOUBLE]
        if not double_bonds:
            return False, "Fails structural check: Not an alkene."

        five_membered_rings_atoms = [set(ring) for ring in ring_info.AtomRings() if len(ring) == 5]
        
        surviving_cyclopentane_found = False
        for ring_atoms in five_membered_rings_atoms:
            # A ring survives if its internal bonds are all single bonds.
            # A double bond is "in" a ring if both its atoms are members of that ring.
            is_db_in_this_ring = False
            for db in double_bonds:
                if db.GetBeginAtomIdx() in ring_atoms and db.GetEndAtomIdx() in ring_atoms:
                    is_db_in_this_ring = True
                    break
            if not is_db_in_this_ring:
                surviving_cyclopentane_found = True
                break
        
        if not surviving_cyclopentane_found:
            return False, "Fails mechanism check: The double bond is inside the 5-membered ring, so it would be opened, not preserved."

        # If all checks pass, the candidate is valid.
        return True, "Satisfies all constraints for the ROCM reaction."

    # --- Main Logic ---
    # 1. Check if the proposed answer is a valid candidate.
    proposed_candidate = options.get(llm_answer_key)
    if not proposed_candidate:
        return f"Incorrect. The answer key '{llm_answer_key}' is not a valid option."

    is_valid, reason = check_candidate(proposed_candidate["smiles"])

    if not is_valid:
        return f"Incorrect. The proposed answer {llm_answer_key} ({proposed_candidate['name']}) is wrong because: {reason}"

    # 2. Verify that all other options fail the checks.
    for key, data in options.items():
        if key == llm_answer_key:
            continue
        passes, reason_for_failure = check_candidate(data["smiles"])
        if passes:
            # This indicates a flaw in the logic or multiple correct answers.
            return f"Incorrect. The logic is flawed because option {key} ({data['name']}) also satisfies all constraints."

    # 3. If the proposed answer is valid and all others are invalid, the answer is correct.
    return "Correct"

# Run the check and print the result.
if __name__ == '__main__':
    # Ensure rdkit is available before running the main logic
    if "rdkit" not in globals() or Chem.MolFromSmiles is None:
        # The error message would have been printed at the import stage.
        pass
    else:
        result = check_correctness()
        print(result)