import sys
import numpy as np
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

def check_cis_esters(mol: Chem.Mol) -> bool | None:
    """
    Checks if the two ester groups are cis to each other by calculating the
    dihedral angle between the hydrogens on their respective attachment carbons.
    A small dihedral angle (< 90 degrees) indicates a cis relationship.
    """
    # Find the two adjacent carbons bearing the ester groups
    patt = Chem.MolFromSmarts('[CX4H1](C(=O)O[CH3])-[CX4H1](C(=O)O[CH3])')
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None  # Pattern not found

    # Get atom indices from the first match
    c1_idx, c2_idx = matches[0][0], matches[0][1]
    c1 = mol.GetAtomWithIdx(c1_idx)
    c2 = mol.GetAtomWithIdx(c2_idx)

    # Find the hydrogens attached to these carbons
    h1_idx, h2_idx = -1, -1
    for neighbor in c1.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            h1_idx = neighbor.GetIdx()
            break
    for neighbor in c2.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:
            h2_idx = neighbor.GetIdx()
            break
    
    if h1_idx == -1 or h2_idx == -1:
        return None  # Hydrogens not found

    conf = mol.GetConformer()
    dihedral = AllChem.GetDihedralDeg(conf, h1_idx, c1_idx, c2_idx, h2_idx)
    
    return abs(dihedral) < 90

def check_anti_isomer(mol: Chem.Mol) -> bool | None:
    """
    Checks if the molecule is the anti-isomer by comparing the relative
    positions of the methano bridge and the ester groups. A negative dot
    product of vectors from the center indicates they are on opposite sides (anti).
    """
    conf = mol.GetConformer()

    # Find the four central cyclobutane atoms (part of 3 rings)
    patt_central = Chem.MolFromSmarts('[c,C;R3]')
    central_atoms_idx = [match[0] for match in mol.GetSubstructMatches(patt_central)]
    if len(central_atoms_idx) != 4:
        return None

    # Find the methano bridge apex atom (CH2 in a 5-membered ring)
    patt_methano = Chem.MolFromSmarts('[CH2;r5]')
    methano_apex_idx_list = [match[0] for match in mol.GetSubstructMatches(patt_methano)]
    if not methano_apex_idx_list:
        return None
    methano_apex_idx = methano_apex_idx_list[0]

    # Find the two ester groups' carbonyl carbons
    patt_esters = Chem.MolFromSmarts('[C;!R](=O)O[CH3]')
    ester_c_idx = [match[0] for match in mol.GetSubstructMatches(patt_esters)]
    if len(ester_c_idx) != 2:
        return None

    # Calculate average positions
    central_pos = np.mean([np.array(conf.GetAtomPosition(i)) for i in central_atoms_idx], axis=0)
    methano_pos = np.array(conf.GetAtomPosition(methano_apex_idx))
    ester_pos = np.mean([np.array(conf.GetAtomPosition(i)) for i in ester_c_idx], axis=0)

    # Create vectors from the center of the molecule
    vec_methano = methano_pos - central_pos
    vec_ester = ester_pos - central_pos

    # Negative dot product means they are on opposite sides (anti)
    return np.dot(vec_methano, vec_ester) < 0

def check_correctness():
    """
    Main function to check the correctness of the provided answer.
    """
    smiles_options = {
        'A': "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'B': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O",
        'C': "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O",
        'D': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O"
    }
    llm_answer = 'A'

    results = {}
    for option, smiles in smiles_options.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            results[option] = {'error': f'Invalid SMILES for option {option}.'}
            continue
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
            results[option] = {'error': f'Could not generate 3D coordinates for option {option}.'}
            continue
        
        is_cis = check_cis_esters(mol)
        is_anti = check_anti_isomer(mol)
        results[option] = {'is_cis': is_cis, 'is_anti': is_anti}

    # The major product must be cis and anti
    predicted_correct_option = None
    for option, res in results.items():
        if res.get('is_cis') and res.get('is_anti'):
            predicted_correct_option = option
            break
    
    if llm_answer == predicted_correct_option:
        return "Correct"
    else:
        llm_res = results.get(llm_answer, {})
        if 'error' in llm_res:
            return f"Incorrect. The SMILES for the provided answer {llm_answer} could not be processed: {llm_res['error']}"
        
        if not llm_res.get('is_cis'):
            return f"Incorrect. The provided answer {llm_answer} is wrong because the major product must have cis-dicarboxylate groups (originating from maleic anhydride), but structure {llm_answer} has trans-dicarboxylate groups."
        
        if not llm_res.get('is_anti'):
            return f"Incorrect. The provided answer {llm_answer} is wrong because the major product is the anti-isomer, formed from the sterically favored attack of cyclopentadiene. Structure {llm_answer} is the syn-isomer, which would be a minor product."
            
        if predicted_correct_option is None:
            return "Incorrect. The provided answer is {llm_answer}, but none of the options satisfy the required chemical principles (cis-diester and anti-isomer)."
        else:
            return f"Incorrect. The provided answer is {llm_answer}, but the major product that satisfies the stereochemical constraints (cis-diester and anti-isomer) is option {predicted_correct_option}."

# Run the check
result = check_correctness()
print(result)