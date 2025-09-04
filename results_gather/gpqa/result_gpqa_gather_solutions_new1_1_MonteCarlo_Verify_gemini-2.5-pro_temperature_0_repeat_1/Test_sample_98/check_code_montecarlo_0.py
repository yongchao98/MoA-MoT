import sys
from io import StringIO

def check_answer():
    """
    Checks the correctness of the LLM's answer by:
    1. Defining the four candidate molecules from the question using SMILES strings.
    2. Defining the required NMR splitting patterns (dtq and dtt).
    3. Using the RDKit library to analyze the neighborhood of each proton in each molecule.
    4. Identifying the single molecule that satisfies both NMR constraints.
    5. Comparing the correct molecule's label with the LLM's provided answer.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Execution Error: RDKit is not installed. Please install it using 'pip install rdkit'."

    # Define the four candidate structures based on the original question's labels.
    # SMILES strings provide an unambiguous representation of the molecules.
    structures = {
        'A': 'CCC(C)C(C)C(=O)O',        # CH3CH2C(H)(CH3)C(H)(CH3)COOH (2,3-dimethylpentanoic acid)
        'B': 'CCC(CC)C(CC)C(=O)O',      # CH3CH2C(H)(C2H5)C(H)(C2H5)COOH (2,3-diethylpentanoic acid)
        'C': 'CC(C)C(C)CC(=O)O',        # CH3C(H)(CH3)C(H)(CH3)CH2COOH (3,4-dimethylpentanoic acid)
        'D': 'CC(CC)C(CC)CC(=O)O'       # CH3C(H)(C2H5)C(H)(C2H5)CH2COOH (3,4-diethylpentanoic acid)
    }

    # Define the required neighbor hydrogen counts for the complex splitting patterns.
    # dtq: doublet (1H neighbor), triplet (2H neighbors), quartet (3H neighbors)
    dtq_neighbors = sorted([1, 2, 3])
    # dtt: doublet (1H neighbor), triplet (2H neighbors), triplet (another 2H neighbors)
    dtt_neighbors = sorted([1, 2, 2])

    analysis_results = {}

    for label, smiles in structures.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        has_dtq_proton = False
        has_dtt_proton = False

        # We analyze the methine (CH) protons, as they are most likely to show complex splitting.
        methine_carbons = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 1 and not atom.IsInRing()]

        for atom in methine_carbons:
            neighbor_h_counts = []
            # Get the number of hydrogens on each adjacent carbon atom.
            # Coupling to the COOH proton is usually not observed or is very broad, so we ignore non-carbon neighbors.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    neighbor_h_counts.append(neighbor.GetTotalNumHs())
            
            # Check if this proton's neighbors match the required patterns
            if sorted(neighbor_h_counts) == dtq_neighbors:
                has_dtq_proton = True
            if sorted(neighbor_h_counts) == dtt_neighbors:
                has_dtt_proton = True
        
        analysis_results[label] = {'has_dtq': has_dtq_proton, 'has_dtt': has_dtt_proton}

    # Identify the molecule that satisfies BOTH conditions
    correct_label = None
    for label, res in analysis_results.items():
        if res['has_dtq'] and res['has_dtt']:
            correct_label = label
            break
    
    # The final answer provided by the LLM being checked
    llm_answer_label = 'A'

    # Final verification
    if correct_label is None:
        return "Analysis Error: No single molecule from the options satisfies both NMR conditions (dtq and dtt). The problem statement or options may be flawed."

    if llm_answer_label == correct_label:
        return "Correct"
    else:
        # The LLM's reasoning identified the correct structure but mislabeled it.
        # The reasoning text says: "A) CH₃C(H)(C₂H₅)C(H)(C₂H₅)CH₂COOH (3,4-diethylpentanoic acid)"
        # This shows the LLM correctly identified the structure but called it 'A' instead of 'D'.
        return (f"Incorrect. The final answer given is <<<A>>>, but the correct answer is <<<{correct_label}>>>.\n\n"
                f"Reasoning: The problem requires finding a molecule with protons that exhibit both a 'doublet of triplets of quartets' (dtq) and a 'doublet of triplets of triplets' (dtt) signal.\n"
                f"My analysis confirms that only structure D, 3,4-diethylpentanoic acid (`{structures['D']}`), has protons that satisfy both conditions:\n"
                f"- The proton at C3 is adjacent to a CH2, a CH, and another CH2 group, resulting in a 'dtt' signal.\n"
                f"- The proton at C4 is adjacent to a CH, a CH3, and a CH2 group, resulting in a 'dtq' signal.\n\n"
                f"The provided answer's reasoning correctly identifies this same molecular structure but incorrectly calls it 'option A'. According to the original question's list, this structure is option D. Therefore, the final output <<<A>>> is wrong.")

# Run the check
result = check_answer()
print(result)