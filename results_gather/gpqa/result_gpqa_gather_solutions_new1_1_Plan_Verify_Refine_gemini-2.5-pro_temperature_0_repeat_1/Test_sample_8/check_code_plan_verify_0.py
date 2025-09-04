import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

def check_mapping(coords1, coords2, types1, types2, tol=0.15):
    """
    Checks if two sets of coordinates represent the same structure by finding
    a one-to-one mapping between atoms of the same type within a tolerance.
    """
    if len(coords1) != len(coords2):
        return False
    
    # A simple greedy matching algorithm is sufficient for this problem.
    # It tries to find a unique partner for each atom in coords1 from coords2.
    used_indices = set()
    for i in range(len(coords1)):
        found_match = False
        for j in range(len(coords2)):
            if j in used_indices:
                continue
            # Check if atom types match and distance is within tolerance
            if types1[i] == types2[j] and np.linalg.norm(coords1[i] - coords2[j]) < tol:
                used_indices.add(j)
                found_match = True
                break
        if not found_match:
            return False
    return True

def get_point_group(mol, name):
    """
    Analyzes the symmetry of an RDKit molecule and returns its likely point group.
    Focuses on distinguishing between C3, C3v, C3h, and D3h.
    """
    try:
        # 1. Generate and optimize conformer
        if mol.GetNumConformers() == 0:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42  # for reproducibility
            # For the flexible borate, sample multiple conformers to find the ground state
            if "borate" in name:
                AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
                res = AllChem.MMFFOptimizeMoleculeConfs(mol)
                if not res: return "Optimization Failed"
                min_energy_idx = np.argmin([e[1] for e in res])
                conf_id = mol.GetConformers()[min_energy_idx].GetId()
            else:  # For rigid molecules, one conformer is sufficient
                conf_id = AllChem.EmbedMolecule(mol, params)
                if conf_id == -1: return "Embedding Failed"
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
        else:
            conf_id = 0

        conf = mol.GetConformer(conf_id)
        
        # 2. Center the molecule and get atom types
        coords = conf.GetPositions()
        center = np.mean(coords, axis=0)
        centered_coords = coords - center
        atom_types = [atom.GetAtomicNum() for atom in mol.GetAtoms()]

        # 3. Find principal axes of inertia
        _, principal_axes = rdMolTransforms.ComputePrincipalMomentsAndAxes(conf)
        
        # 4. Check for C3 axis
        has_c3, c3_axis = False, None
        for axis_candidate in principal_axes:
            angle = 2 * np.pi / 3
            R = rdMolTransforms.GetRotationMatrix(angle, axis_candidate)
            rotated_coords = np.dot(centered_coords, R.T)
            if check_mapping(rotated_coords, centered_coords, atom_types, atom_types):
                has_c3, c3_axis = True, axis_candidate
                break
        
        if not has_c3:
            return "Low Symmetry (No C3)"

        # 5. Check for horizontal mirror plane (sigma_h)
        # A reflection through a plane with normal n is: v' = v - 2*dot(v,n)*n
        reflection_matrix = np.identity(3) - 2 * np.outer(c3_axis, c3_axis)
        reflected_coords = np.dot(centered_coords, reflection_matrix.T)
        has_sigma_h = check_mapping(reflected_coords, centered_coords, atom_types, atom_types)

        # 6. Check for perpendicular C2 axes (to distinguish D3h)
        has_perp_c2 = False
        if has_sigma_h:
            # Find a vector perpendicular to c3_axis to test for a C2 axis
            v_perp = np.array([1.0, 0.0, 0.0])
            if np.abs(np.dot(v_perp, c3_axis)) > 0.95: v_perp = np.array([0.0, 1.0, 0.0])
            c2_candidate_axis = np.cross(c3_axis, v_perp)
            c2_candidate_axis /= np.linalg.norm(c2_candidate_axis)
            
            R_c2 = rdMolTransforms.GetRotationMatrix(np.pi, c2_candidate_axis)
            rotated_coords_c2 = np.dot(centered_coords, R_c2.T)
            if check_mapping(rotated_coords_c2, centered_coords, atom_types, atom_types):
                has_perp_c2 = True

        # 7. Check for vertical mirror planes (to distinguish C3v)
        has_sigma_v = False
        if not has_sigma_h:
            # A vertical plane contains the C3 axis. Find an atom not on the axis to define the plane.
            for atom_pos in centered_coords:
                if np.linalg.norm(np.cross(atom_pos, c3_axis)) > 0.1: # if atom is not on axis
                    plane_normal = np.cross(c3_axis, atom_pos)
                    plane_normal /= np.linalg.norm(plane_normal)
                    
                    reflection_matrix_v = np.identity(3) - 2 * np.outer(plane_normal, plane_normal)
                    reflected_coords_v = np.dot(centered_coords, reflection_matrix_v.T)
                    if check_mapping(reflected_coords_v, centered_coords, atom_types, atom_types):
                        has_sigma_v = True
                        break
        
        # 8. Conclude point group based on found elements
        if has_c3 and has_sigma_h and has_perp_c2: return "D3h"
        if has_c3 and has_sigma_h and not has_perp_c2: return "C3h"
        if has_c3 and has_sigma_v and not has_sigma_h: return "C3v"
        if has_c3: return "C3"
        
        return "Symmetry analysis inconclusive"

    except Exception as e:
        return f"Error during analysis: {e}"

def check_correctness():
    """
    Main function to check the correctness of the proposed answer.
    """
    # The provided answer is B.
    proposed_answer_key = "B"
    
    molecules = {
        "A": ("triisopropyl borate", "CC(C)OB(OC(C)C)OC(C)C"),
        "B": ("triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone", "C1=C2C(=O)OC(=O)C2=C3C4=C(C=C5C(=O)OC(=O)C5=C4)C6=C3C=C1C(=O)OC6=O"),
        "C": ("benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone", "C12=C(C3=C(C4=C1C(=O)OC4=O)C(=O)OC3=O)C(=O)OC2=O"),
        "D": ("quinuclidine", "C1CN2CCC1CC2")
    }

    results = {}
    print("Running computational symmetry analysis...")
    for key, (name, smiles) in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results[key] = "SMILES parsing failed"
            continue
        mol = Chem.AddHs(mol)
        results[key] = get_point_group(mol, name)
        print(f"- Molecule {key} ({name}): Found point group -> {results[key]}")

    # Check if the proposed answer is correct based on the analysis
    c3h_molecule_key = None
    for key, pg in results.items():
        if pg == "C3h":
            if c3h_molecule_key is not None:
                return f"Incorrect. The question is ambiguous as multiple molecules were found to have C3h symmetry: {c3h_molecule_key} and {key}."
            c3h_molecule_key = key
    
    if c3h_molecule_key is None:
        return "Incorrect. No molecule was found to have C3h symmetry based on the analysis."

    if c3h_molecule_key == proposed_answer_key:
        # Final verification against known theory
        if results.get("D") == "C3v" and results.get("C") == "D3h":
             return "Correct"
        else:
            return f"Partially correct. While {c3h_molecule_key} was found to be C3h, other molecules were misidentified (e.g., D should be C3v, C should be D3h), suggesting a potential issue in the analysis."
    else:
        return f"Incorrect. The molecule with C3h symmetry was found to be {c3h_molecule_key}, not {proposed_answer_key}."

# Run the check
print(check_correctness())