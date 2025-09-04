import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
import warnings

# Suppress RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')
warnings.filterwarnings("ignore", category=UserWarning)

def find_closest_atom(pos, coords, types, atom_type, tol=0.3):
    """Helper function to find if a transformed atom overlaps with an original atom."""
    for i, (coord, atype) in enumerate(zip(coords, types)):
        if atype == atom_type and np.linalg.norm(pos - coord) < tol:
            return i
    return -1

def analyze_point_group(mol, conf_id=-1, tol=0.3):
    """
    Analyzes the symmetry of an RDKit molecule's conformer and returns its likely point group.
    This function checks for C3, C3h, C3v, and D3h point groups.
    """
    try:
        if mol.GetNumConformers() == 0:
            return "No Conformer"
        
        conf = mol.GetConformer(conf_id)
        coords = conf.GetPositions()
        atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]
        
        # Create a temporary conformer to avoid modifying the original
        temp_conf = Chem.Conformer(conf)
        com = rdMolTransforms.ComputeCentroid(temp_conf)
        rdMolTransforms.Translate(temp_conf, -com)
        coords = temp_conf.GetPositions()

        _, principal_axes = rdMolTransforms.ComputePrincipalMomentsAndAxes(temp_conf)

        # --- Check for C3 axis ---
        has_c3, c3_axis = False, None
        for axis_idx in range(3):
            axis = principal_axes[:, axis_idx]
            theta = 2 * np.pi / 3
            rotation_matrix = rdMolTransforms.GetRotationMatrix(theta, axis)
            rotated_coords = (coords @ rotation_matrix.T)
            
            mapping = [find_closest_atom(rotated_coords[i], coords, atom_types, atom_types[i], tol) for i in range(len(coords))]
            if -1 not in mapping:
                has_c3, c3_axis = True, axis
                break
        
        if not has_c3: return "Low symmetry (No C3)"

        # --- Check for horizontal mirror plane (sigma_h) ---
        has_sigma_h = False
        if c3_axis is not None:
            reflected_coords = coords - 2 * np.outer(np.dot(coords, c3_axis), c3_axis)
            mapping = [find_closest_atom(reflected_coords[i], coords, atom_types, atom_types[i], tol) for i in range(len(coords))]
            if -1 not in mapping:
                has_sigma_h = True

        # --- Distinguish between D3h, C3h, C3v, C3 ---
        if has_sigma_h:
            # Could be C3h or D3h. Check for a perpendicular C2.
            perp_vec = np.array([1.0, 0.0, 0.0])
            if np.abs(np.dot(perp_vec, c3_axis)) > 0.95:
                perp_vec = np.array([0.0, 1.0, 0.0])
            c2_axis = np.cross(c3_axis, perp_vec)
            c2_axis /= np.linalg.norm(c2_axis)
            
            theta_c2 = np.pi
            rotation_matrix_c2 = rdMolTransforms.GetRotationMatrix(theta_c2, c2_axis)
            rotated_coords_c2 = (coords @ rotation_matrix_c2.T)
            
            mapping_c2 = [find_closest_atom(rotated_coords_c2[i], coords, atom_types, atom_types[i], tol) for i in range(len(coords))]
            if -1 not in mapping_c2:
                return "D3h"
            else:
                return "C3h"
        else:
            # No sigma_h. Could be C3v or C3. Check for a vertical mirror plane (sigma_v).
            for i in range(len(coords)):
                atom_pos = coords[i]
                if np.linalg.norm(np.cross(atom_pos, c3_axis)) < 1e-2: continue
                
                plane_normal = np.cross(c3_axis, atom_pos)
                if np.linalg.norm(plane_normal) < 1e-2: continue
                plane_normal /= np.linalg.norm(plane_normal)
                
                reflected_coords = coords - 2 * np.outer(np.dot(coords, plane_normal), plane_normal)
                mapping = [find_closest_atom(reflected_coords[i], coords, atom_types, atom_types[i], tol) for i in range(len(coords))]
                if -1 not in mapping:
                    return "C3v"
            return "C3"

    except Exception:
        return "Analysis Error"

def check_correctness():
    """
    Checks the correctness of the answer by determining the point group of each molecule.
    """
    smiles_dict = {
        "A": ("benzo[...]trifuran[...]hexaone", "O=C1OC(=O)c2c3c(c4c(c2)C(=O)OC4=O)C(=O)OC3=O"),
        "B": ("triisopropyl borate", "CC(C)OB(OC(C)C)OC(C)C"),
        "C": ("quinuclidine", "C1CN2CCC1CC2"),
        "D": ("triphenyleno[...]trifuran[...]hexaone", "O=C1OC(=O)c2cc3c4cc5c(cc4c(cc21)C(=O)OC3=O)C(=O)OC5=O")
    }
    
    results = {}
    
    for key, (name, smi) in smiles_dict.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results[key] = f"SMILES Parse Error"
            continue
        mol = Chem.AddHs(mol)
        
        if key == "B": # Flexible molecule
            num_confs = 100
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=AllChem.ETKDGv3())
            if not cids:
                results[key] = "Conformer generation failed"
                continue
            
            AllChem.MMFFOptimizeMoleculeConfs(mol)
            
            point_groups_found = {analyze_point_group(mol, conf_id=cid) for cid in cids}
            
            # Determine the highest symmetry point group found
            if "D3h" in point_groups_found: results[key] = "D3h"
            elif "C3v" in point_groups_found: results[key] = "C3v"
            elif "C3h" in point_groups_found: results[key] = "C3h"
            elif "C3" in point_groups_found: results[key] = "C3"
            else: results[key] = "Low symmetry"
        else: # Rigid molecules
            cid = AllChem.EmbedMolecule(mol, params=AllChem.ETKDGv3())
            if cid == -1:
                results[key] = "Conformer generation failed"
                continue
            AllChem.MMFFOptimizeMolecule(mol)
            results[key] = analyze_point_group(mol)

    # --- Verification Logic ---
    # The question asks which molecule *has* C3h symmetry, which is interpreted as
    # belonging to the C3h point group.
    expected_answer = "B"
    c3h_molecule_found = None
    
    for key, pg in results.items():
        if pg == "C3h":
            c3h_molecule_found = key
    
    # Check if the analysis aligns with chemical principles
    if results.get("A") != "D3h":
        return f"Incorrect: Analysis of molecule A is inconsistent. Expected D3h, but code found {results.get('A')}."
    if results.get("C") != "C3v":
        return f"Incorrect: Analysis of molecule C is inconsistent. Expected C3v, but code found {results.get('C')}."
    if results.get("D") != "D3h":
        return f"Incorrect: Analysis of molecule D is inconsistent. Expected D3h, but code found {results.get('D')}."

    # Now, evaluate the final answer
    if c3h_molecule_found == expected_answer:
        return "Correct"
    elif c3h_molecule_found is None:
        return f"Incorrect: The provided answer is {expected_answer}, but no molecule was found to have C3h symmetry. The determined point groups are: A={results.get('A')}, B={results.get('B')}, C={results.get('C')}, D={results.get('D')}."
    else:
        return f"Incorrect: The provided answer is {expected_answer}, but the code identified molecule {c3h_molecule_found} as having C3h symmetry."

# Execute the check and print the result
result = check_correctness()
print(result)