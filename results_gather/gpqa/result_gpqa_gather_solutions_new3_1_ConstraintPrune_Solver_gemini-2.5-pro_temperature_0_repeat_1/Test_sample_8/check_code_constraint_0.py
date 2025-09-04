import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

def find_match(pos, coords, types, atom_type, used_indices, tol=0.25):
    """Helper function to find an unused atom of the correct type close to a given position."""
    for i, (coord, atype) in enumerate(zip(coords, types)):
        if i not in used_indices and atype == atom_type and np.linalg.norm(pos - coord) < tol:
            return i
    return -1

def get_point_group(mol, name):
    """
    Analyzes the symmetry of an RDKit molecule to determine its point group
    from the set {C3, C3v, C3h, D3h, Other}.
    """
    try:
        # 1. Generate and optimize a 3D conformer
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # for reproducibility
        if "triisopropyl borate" in name:
            # For flexible molecules, find the lowest energy conformer
            AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
            res = AllChem.MMFFOptimizeMoleculeConfs(mol)
            if not res: return "Optimization Failed"
            min_energy_idx = np.argmin([e[1] for e in res])
            conf_id = mol.GetConformers()[min_energy_idx].GetId()
        else:
            # For rigid molecules, a single conformer is sufficient
            conf_id = AllChem.EmbedMolecule(mol, params)
            if conf_id == -1: return "Embedding Failed"
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

        conf = mol.GetConformer(conf_id)
        coords = conf.GetPositions()
        center = np.mean(coords, axis=0)
        coords -= center  # Center the molecule
        atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]

        # 2. Find principal axes of inertia as candidates for symmetry axes
        inertia_tensor = np.zeros((3, 3))
        for i, atom in enumerate(mol.GetAtoms()):
            mass = atom.GetMass()
            x, y, z = coords[i]
            inertia_tensor[0, 0] += mass * (y**2 + z**2); inertia_tensor[1, 1] += mass * (x**2 + z**2); inertia_tensor[2, 2] += mass * (x**2 + y**2)
            inertia_tensor[0, 1] -= mass * x * y; inertia_tensor[1, 0] -= mass * x * y
            inertia_tensor[0, 2] -= mass * x * z; inertia_tensor[2, 0] -= mass * x * z
            inertia_tensor[1, 2] -= mass * y * z; inertia_tensor[2, 1] -= mass * y * z
        _, eigvecs = np.linalg.eigh(inertia_tensor)
        
        # 3. Check for a C3 axis
        has_c3, c3_axis = False, None
        for axis in [eigvecs[:, i] for i in range(3)]:
            theta = 2 * np.pi / 3
            K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
            rotated_coords = coords @ R.T
            
            used_indices = set()
            is_match = all((match_idx := find_match(rotated_coords[i], coords, atom_types, atom_types[i], used_indices)) != -1 and not used_indices.add(match_idx) for i in range(len(coords)))
            if is_match:
                has_c3, c3_axis = True, axis
                break
        
        if not has_c3: return "No C3 axis"

        # 4. Check for a horizontal mirror plane (sigma_h)
        n = c3_axis
        reflection_matrix = np.eye(3) - 2 * np.outer(n, n)
        reflected_coords = coords @ reflection_matrix.T
        used_indices = set()
        has_sigma_h = all((match_idx := find_match(reflected_coords[i], coords, atom_types, atom_types[i], used_indices)) != -1 and not used_indices.add(match_idx) for i in range(len(coords)))

        if has_sigma_h:
            # 5. Check for a perpendicular C2 axis (to distinguish from D3h)
            perp_vec = np.array([1.0, 0.0, 0.0])
            if np.abs(np.dot(perp_vec, c3_axis)) > 0.95: perp_vec = np.array([0.0, 1.0, 0.0])
            c2_axis = np.cross(c3_axis, perp_vec); c2_axis /= np.linalg.norm(c2_axis)
            
            theta_c2 = np.pi
            K_c2 = np.array([[0, -c2_axis[2], c2_axis[1]], [c2_axis[2], 0, -c2_axis[0]], [-c2_axis[1], c2_axis[0], 0]])
            R_c2 = np.eye(3) + np.sin(theta_c2) * K_c2 + (1 - np.cos(theta_c2)) * (K_c2 @ K_c2)
            rotated_coords_c2 = coords @ R_c2.T
            used_indices = set()
            has_perp_c2 = all((match_idx := find_match(rotated_coords_c2[i], coords, atom_types, atom_types[i], used_indices)) != -1 and not used_indices.add(match_idx) for i in range(len(coords)))
            
            return "D3h" if has_perp_c2 else "C3h"
        else:
            # 6. Check for a vertical mirror plane (sigma_v) for C3v
            atom_for_plane = coords[np.argmax(np.linalg.norm(np.cross(coords, c3_axis), axis=1))]
            plane_normal = np.cross(c3_axis, atom_for_plane)
            if np.linalg.norm(plane_normal) > 1e-4:
                plane_normal /= np.linalg.norm(plane_normal)
                reflection_matrix_v = np.eye(3) - 2 * np.outer(plane_normal, plane_normal)
                reflected_coords_v = coords @ reflection_matrix_v.T
                used_indices = set()
                has_sigma_v = all((match_idx := find_match(reflected_coords_v[i], coords, atom_types, atom_types[i], used_indices)) != -1 and not used_indices.add(match_idx) for i in range(len(coords)))
                if has_sigma_v: return "C3v"
            return "C3"
    except Exception:
        return "Analysis Error"

def check_correctness():
    """Main function to run the analysis and check the answer."""
    # The question's options are:
    # A) quinuclidine
    # B) triphenyleno[...]
    # C) triisopropyl borate
    # D) benzo[...]
    molecules = {
        "A) quinuclidine": "C1CN2CCC1CC2",
        "B) triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone": "C1=C2C(=O)OC(=O)C2=C3C4=C(C=C5C(=O)OC(=O)C5=C4)C6=C3C=C1C(=O)OC6=O",
        "C) triisopropyl borate": "CC(C)OB(OC(C)C)OC(C)C",
        "D) benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone": "C12=C(C3=C(C4=C1C(=O)OC4=O)C(=O)OC3=O)C(=O)OC2=O",
    }
    
    results = {}
    print("--- Computational Symmetry Analysis ---")
    for name, smi in molecules.items():
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        point_group = get_point_group(mol, name)
        results[name] = point_group
        print(f"Result for {name}: {point_group}")
    print("------------------------------------")

    # The provided answer is 'B'. Let's verify.
    correct_answer_key = "B) triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone"
    
    # Check 1: Does the chosen answer have C3h symmetry?
    if results.get(correct_answer_key) != "C3h":
        return f"Incorrect. The provided answer (B) was computationally determined to have point group '{results.get(correct_answer_key)}', not C3h."

    # Check 2: Is it the *only* molecule with C3h symmetry?
    for name, pg in results.items():
        if name != correct_answer_key and pg == "C3h":
            return f"Incorrect. Another molecule ({name}) was also found to have C3h symmetry, making the answer ambiguous."

    # Check 3: Do the other molecules have the expected point groups, validating the method?
    # A molecule with D3h symmetry also contains C3h elements, but its full point group is D3h.
    # In such questions, the answer is the molecule whose highest point group is the one specified.
    if results.get("D) benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone") == "D3h" and \
       results.get("A) quinuclidine") == "C3v" and \
       results.get("C) triisopropyl borate") == "C3":
        # All checks pass. The chosen answer is the only one with exactly C3h symmetry.
        return "Correct"
    else:
        return f"The answer appears correct, but the analysis of other molecules produced unexpected results, which may indicate a method limitation. Expected A=C3v, C=C3, D=D3h. Found A={results.get('A) quinuclidine')}, C={results.get('C) triisopropyl borate')}, D={results.get('D) benzo[1,2-c:3,4-c'']trifuran-1,3,4,6,7,9-hexaone')}."

# Run the check
result = check_correctness()
print(result)