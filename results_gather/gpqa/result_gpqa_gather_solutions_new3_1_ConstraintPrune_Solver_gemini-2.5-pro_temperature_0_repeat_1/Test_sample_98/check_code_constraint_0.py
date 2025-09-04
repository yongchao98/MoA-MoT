import sys
from collections import Counter

# This script requires the RDKit library.
# You can install it via pip: pip install rdkit
# Or via conda: conda install -c conda-forge rdkit

try:
    from rdkit import Chem
except ImportError:
    print("RDKit library not found. Please install it to run this check:")
    print("pip install rdkit")
    sys.exit(1)

def check_nmr_splitting_patterns():
    """
    Analyzes four candidate molecules to determine which ones match the
    ¹H NMR splitting patterns described in the problem.
    - dtq: doublet of triplets of quartets
    - dtt: doublet of triplets of triplets
    """
    # Define the structures from the question using SMILES strings
    # A) CH3C(H)(C2H5)C(H)(C2H5)CH2COOH -> 3,4-diethylpentanoic acid
    #    Structure: COOH-CH2-CH(Et)-CH(Et)-CH3
    smiles_A = "CCC(CC)C(CC)CC(=O)O"
    
    # B) CH3C(H)(CH3)C(H)(CH3)CH2COOH -> 3,4-dimethylpentanoic acid
    #    Structure: COOH-CH2-CH(Me)-CH(Me)-CH3
    smiles_B = "CC(C)C(C)CC(=O)O"
    
    # C) CH3CH2C(H)(C2H5)C(H)(C2H5)COOH -> 2,3-diethylpentanoic acid
    #    Structure: COOH-CH(Et)-CH(Et)-CH2-CH3
    smiles_C = "CCC(CC)C(CC)C(=O)O"
    
    # D) CH3CH2C(H)(CH3)C(H)(CH3)COOH -> 2,3-dimethylpentanoic acid
    #    Structure: COOH-CH(Me)-CH(Me)-CH2-CH3
    smiles_D = "CCC(C)C(C)C(=O)O"

    structures = {
        "A": smiles_A,
        "B": smiles_B,
        "C": smiles_C,
        "D": smiles_D,
    }

    # Define the required splitting patterns based on neighbor H counts
    # Using Counter to be order-independent
    dtq_pattern = Counter([1, 2, 3])
    dtt_pattern = Counter([1, 2, 2])

    analysis_results = {}

    for label, smiles in structures.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            analysis_results[label] = "Invalid SMILES string"
            continue
        
        mol = Chem.AddHs(mol) # Add explicit hydrogens to the graph

        has_dtq = False
        has_dtt = False

        # Iterate through each atom to find methine (CH) groups,
        # as these are the most likely to have complex splitting.
        for atom in mol.GetAtoms():
            # We are interested in a proton on a carbon atom (a methine group)
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 1:
                
                neighbor_h_counts = []
                # Get its neighboring atoms to check for coupling
                for neighbor in atom.GetNeighbors():
                    # We only consider coupling to protons on adjacent carbons
                    if neighbor.GetSymbol() == 'C':
                        # Count the hydrogens on this neighboring carbon
                        neighbor_h_counts.append(neighbor.GetTotalNumHs())
                
                # Check if the neighbor counts match our patterns
                if Counter(neighbor_h_counts) == dtq_pattern:
                    has_dtq = True
                if Counter(neighbor_h_counts) == dtt_pattern:
                    has_dtt = True
        
        analysis_results[label] = {"has_dtq": has_dtq, "has_dtt": has_dtt}

    # Determine which molecule(s) satisfy the conditions
    correct_molecules = []
    for label, res in analysis_results.items():
        if isinstance(res, dict) and res["has_dtq"] and res["has_dtt"]:
            correct_molecules.append(label)
            
    # The final answer provided in the prompt is 'A'
    provided_answer = "A"

    if not correct_molecules:
        return f"Incorrect. The provided answer is {provided_answer}, but the analysis found no molecule that satisfies both conditions (having both a 'dtq' and a 'dtt' signal). There might be an error in the problem statement or the options. Analysis details: {analysis_results}"
    
    if provided_answer not in correct_molecules:
        return f"Incorrect. The provided answer is {provided_answer}, but the analysis shows that molecule(s) {correct_molecules} are the ones that satisfy the conditions. Molecule {provided_answer} does not have both required signals. Analysis details: {analysis_results}"

    if len(correct_molecules) > 1:
         return f"Ambiguous. The provided answer {provided_answer} is one of the correct options, but the analysis found multiple molecules ({correct_molecules}) that satisfy the conditions. Analysis details: {analysis_results}"

    return "Correct"

# Run the check
result = check_nmr_splitting_patterns()
print(result)