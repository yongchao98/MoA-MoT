import sys
from io import StringIO

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the chemistry question.
    It uses the RDKit library to parse chemical structures and analyze proton environments
    to verify the predicted 1H NMR splitting patterns.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library not found. Please install it using 'pip install rdkit'."

    # --- Problem Definition ---
    # The question asks to identify a compound that has both:
    # 1. A proton signal that is a doublet of triplets of quartets (dtq).
    #    This means a CH proton is adjacent to a CH, a CH2, and a CH3 group.
    #    Neighbor H counts: {1, 2, 3}
    # 2. A proton signal that is a doublet of triplets of triplets (dtt).
    #    This means a CH proton is adjacent to a CH and two different CH2 groups.
    #    Neighbor H counts: {1, 2, 2}

    # Condensed formulas from the question and their corresponding SMILES strings
    molecules = {
        "A": ("CH3CH2C(H)(CH3)C(H)(CH3)COOH", "CCC(C)C(C)C(=O)O"),      # 2,3-dimethylpentanoic acid
        "B": ("CH3CH2C(H)(C2H5)C(H)(C2H5)COOH", "CCC(CC)C(CC)C(=O)O"),    # 2,3-diethylpentanoic acid
        "C": ("CH3C(H)(CH3)C(H)(CH3)CH2COOH", "CC(C)C(C)CC(=O)O"),      # 3,4-dimethylpentanoic acid
        "D": ("CH3C(H)(C2H5)C(H)(C2H5)CH2COOH", "CC(CC)C(CC)CC(=O)O")     # 3,4-diethylpentanoic acid
    }
    
    llm_answer = "D"

    def get_neighbor_h_counts(mol, atom):
        """
        For a given C-H proton (represented by its carbon atom), find the number
        of protons on each adjacent carbon atom.
        """
        # We are only interested in methine (CH) protons for these complex signals
        if not (atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 1):
            return None

        neighbor_h_counts = []
        for neighbor in atom.GetNeighbors():
            # We only care about coupling to protons on adjacent carbons
            if neighbor.GetSymbol() == 'C':
                num_hs = neighbor.GetTotalNumHs()
                if num_hs > 0:
                    neighbor_h_counts.append(num_hs)
        
        # We need exactly 3 different neighboring groups of protons for dtq or dtt
        if len(neighbor_h_counts) != 3:
            return None
            
        return sorted(neighbor_h_counts)

    def analyze_molecule(smiles):
        """
        Checks a molecule for the presence of protons that would give
        a dtq or dtt signal.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False, False
        mol = Chem.AddHs(mol)

        has_dtq = False
        has_dtt = False

        for atom in mol.GetAtoms():
            neighbor_counts = get_neighbor_h_counts(mol, atom)
            
            if neighbor_counts is None:
                continue

            # Check for dtq: neighbors are 1H, 2H, 3H
            if neighbor_counts == [1, 2, 3]:
                has_dtq = True
            
            # Check for dtt: neighbors are 1H, 2H, 2H
            if neighbor_counts == [1, 2, 2]:
                has_dtt = True
                
        return has_dtq, has_dtt

    # --- Verification Logic ---
    analysis_results = {}
    for label, (name, smiles) in molecules.items():
        has_dtq, has_dtt = analyze_molecule(smiles)
        analysis_results[label] = {"dtq": has_dtq, "dtt": has_dtt}

    correct_candidates = []
    for label, res in analysis_results.items():
        if res["dtq"] and res["dtt"]:
            correct_candidates.append(label)

    # --- Final Judgement ---
    if len(correct_candidates) == 0:
        return "Incorrect. The analysis shows that no single molecule from the options satisfies both NMR conditions (having a dtq and a dtt signal). There might be an error in the question's premise or the provided options."
    
    if len(correct_candidates) > 1:
        return f"Incorrect. The question is ambiguous as molecules {', '.join(correct_candidates)} both satisfy the NMR conditions."

    # Exactly one correct candidate was found
    true_answer = correct_candidates[0]
    if llm_answer == true_answer:
        return "Correct"
    else:
        llm_res = analysis_results[llm_answer]
        reason = f"Incorrect. The provided answer is {llm_answer}, but the only molecule that satisfies both NMR conditions is {true_answer} ({molecules[true_answer][0]}).\n"
        reason += f"Analysis of molecule {llm_answer} shows it has a dtq signal: {llm_res['dtq']} and a dtt signal: {llm_res['dtt']}. It does not meet both criteria."
        return reason

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

# Run the checking function
result = check_correctness()

# Restore stdout
sys.stdout = old_stdout

# Print the final result
print(result)