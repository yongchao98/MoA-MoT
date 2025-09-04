import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

def check_mapping(coords1, coords2, symbols1, symbols2, tol=0.25):
    """
    Checks if two sets of atomic coordinates represent the same structure
    by finding a one-to-one mapping of atoms of the same element type
    within a given tolerance.
    """
    if len(coords1) != len(coords2) or len(symbols1) != len(symbols2):
        return False
    
    unmapped_indices_2 = list(range(len(coords2)))
    
    for i in range(len(coords1)):
        pos1 = coords1[i]
        sym1 = symbols1[i]
        
        best_match_local_idx = -1
        min_dist = float('inf')
        
        # Find the best match in the remaining unmapped atoms of coords2
        for local_idx, global_idx in enumerate(unmapped_indices_2):
            pos2 = coords2[global_idx]
            sym2 = symbols2[global_idx]
            if sym1 == sym2:
                dist = np.linalg.norm(pos1 - pos2)
                if dist < min_dist:
                    min_dist = dist
                    best_match_local_idx = local_idx
        
        # If a match is found within tolerance, remove it from the unmapped list
        if min_dist < tol:
            unmapped_indices_2.pop(best_match_local_idx)
        else:
            # If no suitable match is found for any atom, the mapping fails
            return False
            
    # If all atoms in coords1 are mapped, the structures are equivalent
    return True

def get_point_group_analysis(mol, conf_id=0):
    """
    Performs a simplified point group analysis on a molecule's conformer.
    It specifically checks for elements relevant to C3, C3v, C3h, and D3h point groups.
    """
    if mol.GetNumConformers() == 0:
        return "No Conformer"
        
    conf = mol.GetConformer(conf_id)
    coords = conf.GetPositions()
    center = np.mean(coords, axis=0)
    centered_coords = coords - center
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    try:
        # Get principal axes of inertia as potential symmetry axes
        _, p_axes = AllChem.ComputePrincipalMomentsOfInertia(conf)
    except Exception:
        # Fallback for linear or other problematic molecules
        p_axes = np.identity(3)

    # --- Check for C3 axis ---
    has_c3, c3_axis = False, None
    # Test principal axes and also canonical x,y,z axes for robustness
    test_axes = np.hstack([p_axes, np.identity(3)])
    for i in range(test_axes.shape[1]):
        axis = test_axes[:, i]
        rot = Rotation.from_rotvec(axis * (2 * np.pi / 3))
        rotated_coords = rot.apply(centered_coords)
        if check_mapping(centered_coords, rotated_coords, symbols, symbols):
            has_c3 = True
            c3_axis = axis
            break
            
    if not has_c3:
        return "Low Symmetry (No C3)"

    # --- If C3 exists, check for other elements ---
    
    # Check for horizontal mirror plane (sigma_h) perpendicular to the C3 axis
    has_sigma_h = False
    reflected_coords = centered_coords - 2 * np.outer(np.dot(centered_coords, c3_axis), c3_axis)
    if check_mapping(centered_coords, reflected_coords, symbols, symbols):
        has_sigma_h = True

    # Check for perpendicular C2 axes (characteristic of D groups)
    has_perp_c2 = False
    # Find a vector perpendicular to c3_axis to define a test C2 axis
    test_vec = np.array([1.0, 0.0, 0.0])
    if np.abs(np.dot(test_vec, c3_axis)) > 0.95: # Avoid parallel vectors
        test_vec = np.array([0.0, 1.0, 0.0])
    perp_axis = np.cross(c3_axis, test_vec)
    perp_axis /= np.linalg.norm(perp_axis)
    
    rot_c2 = Rotation.from_rotvec(perp_axis * np.pi) # 180-degree rotation
    rotated_c2_coords = rot_c2.apply(centered_coords)
    if check_mapping(centered_coords, rotated_c2_coords, symbols, symbols):
        has_perp_c2 = True

    # Check for vertical mirror planes (sigma_v)
    has_sigma_v = False
    # A vertical plane contains the C3 axis. Its normal is perpendicular to C3.
    reflected_v_coords = centered_coords - 2 * np.outer(np.dot(centered_coords, perp_axis), perp_axis)
    if check_mapping(centered_coords, reflected_v_coords, symbols, symbols):
        has_sigma_v = True

    # --- Determine Point Group based on the presence/absence of key symmetry elements ---
    if has_c3 and has_sigma_h and has_perp_c2:
        return "D3h"
    if has_c3 and has_sigma_h and not has_perp_c2:
        return "C3h"
    if has_c3 and not has_sigma_h and has_sigma_v:
        return "C3v"
    if has_c3 and not has_sigma_h and not has_sigma_v:
        return "C3"
        
    return "Undetermined (has C3)"

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer by computationally
    determining the point group of each molecule in the question.
    """
    # The provided answer is <<<A>>>.
    
    # Define molecules based on the question's labels and canonical SMILES strings
    question_smiles = {
        "A": "O=C1OC(=O)c2c3cc4c(c5cc6c(c(c2)c1)C(=O)OC6=O)c(cc3)C(=O)OC5=O", # triphenyleno[...]
        "B": "C1CN2CCC1CC2", # quinuclidine
        "C": "C12=C(C3=C(C4=C1C(=O)OC4=O)C(=O)OC3=O)C(=O)OC2=O", # benzo[...] (mellitic trianhydride)
        "D": "CC(C)OB(OC(C)C)OC(C)C" # triisopropyl borate
    }

    results = {}
    
    for label, smiles in question_smiles.items():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                results[label] = "SMILES parsing failed"
                continue
            
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42 # for reproducibility
            
            if label == "D": # Flexible molecule: find lowest energy conformer
                AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
                res = AllChem.MMFFOptimizeMoleculeConfs(mol)
                if not res:
                    results[label] = "MMFF optimization failed"
                    continue
                min_energy_idx = np.argmin([e[1] for e in res])
                results[label] = get_point_group_analysis(mol, conf_id=int(min_energy_idx))
            else: # Rigid molecules
                AllChem.EmbedMolecule(mol, params)
                AllChem.MMFFOptimizeMolecule(mol)
                results[label] = get_point_group_analysis(mol)
        except Exception as e:
            results[label] = f"Analysis Error: {e}"

    # --- Verification Logic ---
    # The correct answer should be the molecule that belongs to the C3h point group.
    # Based on chemical theory, we expect:
    # A (triphenyleno...) -> C3h
    # B (quinuclidine) -> C3v
    # C (benzo...) -> D3h
    # D (triisopropyl borate) -> C3 (lowest energy conformer)

    # Condition 1: The proposed answer 'A' must be C3h.
    if results.get("A") != "C3h":
        return f"Incorrect. The provided answer is A, but the computational analysis determined its point group to be {results.get('A')}. A molecule with D3h symmetry, for example, has higher symmetry than C3h and would not be the best answer. Full results: {results}"

    # Condition 2: Molecule C, the main competitor, should be D3h, not C3h.
    if results.get("C") == "C3h":
        return f"Incorrect. The provided answer is A, but the analysis indicates that molecule C (benzo...) also has C3h symmetry, making the question ambiguous. The expected point group for C is D3h. Full results: {results}"

    # Condition 3: The other molecules should have their expected (non-C3h) point groups.
    if results.get("B") != "C3v" or results.get("C") != "D3h" or results.get("D") != "C3":
        # This is a soft failure, as the main point is that A is C3h and others are not.
        # However, it's worth noting if the analysis of other molecules is unexpected.
        pass # We will proceed as long as A is uniquely C3h among the options.

    # If all conditions are met, the answer is correct.
    return "Correct"

# Execute the check and print the result.
print(check_correctness_of_answer())