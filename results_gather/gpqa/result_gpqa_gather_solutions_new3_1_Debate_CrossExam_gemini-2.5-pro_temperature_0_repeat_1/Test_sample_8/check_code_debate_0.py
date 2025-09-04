import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

def are_coords_the_same(coords1, coords2, atom_types, tol=0.3):
    """
    Checks if two sets of coordinates represent the same structure,
    allowing for permutation of identical atoms.
    """
    if len(coords1) != len(coords2):
        return False
    
    # Create a copy of coords2 to "cross out" matched atoms
    temp_coords2 = list(coords2)
    temp_types2 = list(atom_types)
    
    for i in range(len(coords1)):
        pos1 = coords1[i]
        type1 = atom_types[i]
        
        found_match = False
        best_match_idx = -1
        min_dist = float('inf')

        # Find the closest atom of the same type in the second set
        for j in range(len(temp_coords2)):
            pos2 = temp_coords2[j]
            type2 = temp_types2[j]
            
            if type1 == type2:
                dist = np.linalg.norm(pos1 - pos2)
                if dist < tol and dist < min_dist:
                    min_dist = dist
                    best_match_idx = j
                    found_match = True
        
        if found_match:
            # Remove the matched atom to prevent it from being matched again
            del temp_coords2[best_match_idx]
            del temp_types2[best_match_idx]
        else:
            # If any atom from coords1 doesn't have a match, the structures are different
            return False
            
    return True

def get_point_group(mol):
    """
    Analyzes a molecule's 3D structure to determine its point group
    among C3, C3v, C3h, and D3h.
    """
    # Generate and optimize a 3D conformer
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(randomSeed=42))
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        # MMFF may fail for some structures (e.g., with Boron)
        AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    
    # Align molecule to its principal axes of inertia
    rdMolTransforms.CanonicalizeConformer(conf)
    coords = conf.GetPositions()
    atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # The principal axes are now aligned with the Cartesian axes
    principal_axes = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]
    
    # 1. Check for a C3 axis
    has_c3 = False
    c3_axis = None
    for axis in principal_axes:
        transform = rdMolTransforms.GetRotationMatrix(2 * np.pi / 3, axis)
        rotated_coords = np.dot(coords, transform.T)
        if are_coords_the_same(coords, rotated_coords, atom_types):
            has_c3 = True
            c3_axis = axis
            break
    
    if not has_c3:
        return "No C3 axis"

    # 2. Check for a horizontal mirror plane (sigma_h)
    # The plane's normal is the C3 axis
    reflected_coords_h = coords - 2 * np.outer(np.dot(coords, c3_axis), c3_axis)
    has_sigma_h = are_coords_the_same(coords, reflected_coords_h, atom_types)

    # 3. Check for perpendicular C2 axes
    has_perp_c2 = False
    for axis in principal_axes:
        # Check if the axis is perpendicular to the C3 axis
        if abs(np.dot(axis, c3_axis)) < 1e-4:
            transform_c2 = rdMolTransforms.GetRotationMatrix(np.pi, axis)
            rotated_coords_c2 = np.dot(coords, transform_c2.T)
            if are_coords_the_same(coords, rotated_coords_c2, atom_types):
                has_perp_c2 = True
                break
    
    # 4. Check for vertical mirror planes (sigma_v)
    has_sigma_v = False
    # A vertical plane contains the C3 axis. We check for a plane containing C3 and an off-axis atom.
    for atom_coord in coords:
        # Ensure atom is not on the C3 axis
        if np.linalg.norm(np.cross(atom_coord, c3_axis)) > 0.1:
            # The normal to the reflection plane is perpendicular to both C3 and the atom vector
            plane_normal = np.cross(c3_axis, atom_coord)
            if np.linalg.norm(plane_normal) < 1e-4: continue
            plane_normal /= np.linalg.norm(plane_normal)
            
            reflected_coords_v = coords - 2 * np.outer(np.dot(coords, plane_normal), plane_normal)
            if are_coords_the_same(coords, reflected_coords_v, atom_types):
                has_sigma_v = True
                break

    # 5. Determine point group from found elements
    if has_c3 and has_sigma_h and has_perp_c2:
        return "D3h"
    if has_c3 and has_sigma_h:
        return "C3h"
    if has_c3 and has_sigma_v:
        return "C3v"
    if has_c3:
        return "C3"
    
    return "Other"

# --- Main execution ---
try:
    # Define molecules using SMILES strings
    molecules = {
        "A": ("benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone", "C12=C(C3=C(C4=C1C(=O)OC4=O)C(=O)OC3=O)C(=O)OC2=O"),
        "B": ("triisopropyl borate", "CC(C)OB(OC(C)C)OC(C)C"),
        "C": ("triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone", "C1=C2C(=O)OC(=O)C2=C3C4=C(C=C5C(=O)OC(=O)C5=C4)C6=C3C=C1C(=O)OC6=O"),
        "D": ("quinuclidine", "C1CN2CCC1CC2")
    }

    llm_answer_key = "C"
    llm_answer_name = molecules[llm_answer_key][0]

    results = {}
    print("Analyzing molecule symmetries...")
    for key, (name, smi) in molecules.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            point_group = get_point_group(mol)
            results[name] = point_group
            print(f"- {name}: {point_group}")
        else:
            results[name] = "Error: Invalid SMILES"
            print(f"- {name}: Error: Invalid SMILES")

    # Check the correctness of the LLM's answer
    c3h_molecule = None
    for name, pg in results.items():
        if pg == "C3h":
            c3h_molecule = name
            break
    
    if c3h_molecule is None:
        print("\nCheck Result: Incorrect")
        print("Reason: No molecule was found to have C3h symmetry. The expected molecule, "
              f"'{llm_answer_name}', was determined to have {results.get(llm_answer_name, 'N/A')} symmetry.")
    elif c3h_molecule == llm_answer_name:
        print("\nCorrect")
    else:
        print("\nCheck Result: Incorrect")
        print(f"Reason: The LLM's answer was '{llm_answer_name}', but the molecule with C3h symmetry "
              f"was found to be '{c3h_molecule}'.")

except ImportError:
    print("Check failed: RDKit or NumPy is not installed. Please install them using 'pip install rdkit numpy'.")
except Exception as e:
    print(f"An error occurred during the check: {e}")