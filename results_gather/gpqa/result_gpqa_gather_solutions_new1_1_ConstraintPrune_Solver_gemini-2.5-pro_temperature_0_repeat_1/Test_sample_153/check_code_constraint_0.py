import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def check_answer():
    """
    This function checks the correctness of the identified chemical structure based on spectral data.
    It analyzes four candidate molecules against Mass Spec, IR, and 1H NMR data.
    """
    # Define the candidate molecules using their SMILES strings
    candidates = {
        "A": {"name": "Phenyl chloroformate", "smiles": "O=C(OC1=CC=CC=C1)Cl"},
        "B": {"name": "3-Chloro-2-hydroxybenzaldehyde", "smiles": "C1=CC(=C(C(=C1)C=O)O)Cl"},
        "C": {"name": "2-chlorobenzoic acid", "smiles": "C1=CC=C(C(=C1)C(=O)O)Cl"},
        "D": {"name": "4-chlorobenzoic acid", "smiles": "C1=CC(=CC=C1C(=O)O)Cl"}
    }

    # The final answer provided by the LLM to be checked
    llm_answer_key = "D"

    # --- Define constraints from the spectral data ---
    
    # 1. Mass Spectrometry Constraints
    # M+ peak at m/z = 156, M+2 at 32% -> contains one Cl atom, nominal mass is 156
    def check_ms(mol):
        # Check for exactly one chlorine atom
        formula = rdMolDescriptors.CalcMolFormula(mol)
        if formula.count('Cl') != 1 or 'Cl2' in formula:
            return False, "Does not contain exactly one chlorine atom."
        
        # Calculate nominal mass (integer mass of most abundant isotopes)
        # The M+ peak corresponds to the molecule with 35Cl
        nominal_mass = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'Cl':
                nominal_mass += 35  # Use 35 for the M+ peak isotope
            else:
                nominal_mass += atom.GetMass() # RDKit's GetMass() gives integer mass
        
        if nominal_mass != 156:
            return False, f"Incorrect nominal mass ({nominal_mass}). Expected 156."
        
        return True, "Passes MS check."

    # 2. IR and NMR Functional Group Constraints
    # Broad peak 3500-2700 cm-1 and peak at 1720 cm-1 -> Carboxylic acid
    # NMR peak at 11.0 ppm -> Carboxylic acid proton
    def check_functional_group(mol):
        carboxylic_acid_smarts = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        if mol.HasSubstructMatch(carboxylic_acid_smarts):
            return True, "Contains a carboxylic acid group."
        else:
            return False, "Does not contain a carboxylic acid group, which is required by IR and NMR data."

    # 3. 1H NMR Aromatic Pattern Constraint
    # Two doublets, each integrating to 2H -> 1,4-disubstituted (para) benzene ring
    def check_aromatic_pattern(mol):
        # Pattern for a benzene ring
        benzene_smarts = Chem.MolFromSmarts('c1ccccc1')
        if not mol.HasSubstructMatch(benzene_smarts):
            return False, "Molecule does not contain a benzene ring."

        ring_atoms_indices = mol.GetSubstructMatch(benzene_smarts)
        
        substituent_positions = []
        for i, atom_idx in enumerate(ring_atoms_indices):
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                # If a neighbor is not part of the ring, it's a substituent
                if neighbor.GetIdx() not in ring_atoms_indices:
                    substituent_positions.append(i)
                    break 
        
        # The NMR shows 4 aromatic protons, so the ring must be disubstituted
        if len(substituent_positions) != 2:
            return False, f"Is {len(substituent_positions)}-substituted, but NMR indicates a disubstituted ring (4H)."

        # Check for para (1,4) substitution. Positions are 0-5.
        # Para positions have a difference of 3 (e.g., 0 and 3; 1 and 4; 2 and 5).
        pos_diff = abs(substituent_positions[0] - substituent_positions[1])
        if pos_diff == 3:
            return True, "Is para-substituted, matching the NMR pattern."
        else:
            pattern = "ortho" if pos_diff == 1 or pos_diff == 5 else "meta"
            return False, f"Is {pattern}-substituted, which contradicts the para-pattern (two doublets, 2H each) seen in the NMR."

    # --- Analysis Loop ---
    analysis_results = {}
    for key, data in candidates.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            analysis_results[key] = "Invalid SMILES string."
            continue
        
        mol = Chem.AddHs(mol)

        # Apply checks sequentially
        ms_ok, ms_reason = check_ms(mol)
        if not ms_ok:
            analysis_results[key] = f"Failed MS check: {ms_reason}"
            continue

        fg_ok, fg_reason = check_functional_group(mol)
        if not fg_ok:
            analysis_results[key] = f"Failed functional group check: {fg_reason}"
            continue

        nmr_ok, nmr_reason = check_aromatic_pattern(mol)
        if not nmr_ok:
            analysis_results[key] = f"Failed NMR aromatic check: {nmr_reason}"
            continue
        
        # If all checks pass, this is a valid candidate
        analysis_results[key] = "All spectral data are consistent with this structure."

    # --- Final Verdict ---
    passing_candidates = [k for k, v in analysis_results.items() if "All spectral data" in v]

    if not passing_candidates:
        return "Error: The code found that no candidate matches all the spectral data."
    
    if len(passing_candidates) > 1:
        return f"Error: The code found multiple candidates ({', '.join(passing_candidates)}) that match the data."

    # The single candidate that passed all checks
    correct_key = passing_candidates[0]

    if correct_key == llm_answer_key:
        return "Correct"
    else:
        llm_answer_name = candidates[llm_answer_key]['name']
        correct_answer_name = candidates[correct_key]['name']
        reason_for_failure = analysis_results[llm_answer_key]
        
        return (f"Incorrect. The provided answer is {llm_answer_key} ({llm_answer_name}), but the analysis "
                f"identifies {correct_key} ({correct_answer_name}) as the correct structure.\n"
                f"Reason: The provided answer {llm_answer_key} is incorrect because it {reason_for_failure.lower()}")

# Run the check and print the result
try:
    # rdkit is a required dependency.
    # You can install it via pip: pip install rdkit
    print(check_answer())
except ImportError:
    print("Error: RDKit library not found. Please install it using 'pip install rdkit' to run this check.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
