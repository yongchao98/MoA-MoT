import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def check_organic_synthesis_answer():
    """
    Checks the correctness of the provided answer by verifying the structure against
    the stereochemical constraints of the reaction sequence.
    """
    options = {
        'A': {'smiles': 'O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O'},
        'B': {'smiles': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O'},
        'C': {'smiles': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O'},
        'D': {'smiles': 'O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O'}
    }
    llm_answer = 'C'

    # Helper function to get coordinates
    def get_coords(conf, atom_indices):
        if isinstance(atom_indices, int):
            return np.array(conf.GetAtomPosition(atom_indices))
        coords = [np.array(conf.GetAtomPosition(i)) for i in atom_indices]
        return np.mean(coords, axis=0)

    for key, data in options.items():
        smiles = data['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            data['error'] = "Invalid SMILES string."
            continue
        
        mol = Chem.AddHs(mol)
        # Generate a 3D conformation
        embed_params = AllChem.ETKDGv3()
        embed_params.randomSeed = 42 # for reproducibility
        if AllChem.EmbedMolecule(mol, embed_params) == -1:
            data['error'] = "Failed to generate 3D conformation."
            continue
        AllChem.MMFFOptimizeMolecule(mol)
        AllChem.AssignStereochemistryFrom3D(mol)
        conf = mol.GetConformer()

        # --- Constraint 1: Cis-Diester ---
        patt_ester_carbon = Chem.MolFromSmarts('[CX4H](C(=O)OC)')
        matches = mol.GetSubstructMatches(patt_ester_carbon)
        if len(matches) != 2:
            data['error'] = "Could not find exactly two ester-bearing carbons."
            continue
        
        ester_carbons = [m[0] for m in matches]
        chiral_centers = dict(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        tag1 = chiral_centers.get(ester_carbons[0])
        tag2 = chiral_centers.get(ester_carbons[1])
        data['is_cis'] = (tag1 != tag2)

        # --- Constraint 2: Anti-Addition ---
        sssr = Chem.GetSymmSSSR(mol)
        four_membered_ring = next((r for r in sssr if len(r) == 4), None)
        if not four_membered_ring:
            data['error'] = "Could not find central 4-membered ring."
            continue
        
        patt_methano = Chem.MolFromSmarts('[CH2;r5]') # CH2 in a 5-membered ring
        methano_match = mol.GetSubstructMatches(patt_methano)
        if not methano_match:
            data['error'] = "Could not find methano bridge."
            continue
        methano_idx = methano_match[0][0]

        ring_coords = [get_coords(conf, i) for i in four_membered_ring]
        ring_centroid = np.mean(ring_coords, axis=0)
        normal = np.cross(ring_coords[1] - ring_coords[0], ring_coords[2] - ring_coords[0])
        
        vec_methano = get_coords(conf, methano_idx) - ring_centroid
        vec_ester_bridge = get_coords(conf, ester_carbons) - ring_centroid
        
        data['is_anti'] = (np.dot(vec_methano, normal) * np.dot(vec_ester_bridge, normal) < 0)

        # --- Constraint 3: Endo-Selectivity ---
        patt_cc_double = Chem.MolFromSmarts('[#6]=[#6]')
        norbornene_cc = next((b for b in mol.GetSubstructMatches(patt_cc_double) if mol.GetAtomWithIdx(b[0]).IsInRingSize(5)), None)
        if not norbornene_cc:
            data['error'] = "Could not find norbornene C=C bond."
            continue
        
        dist_to_cc = np.linalg.norm(get_coords(conf, ester_carbons) - get_coords(conf, norbornene_cc))
        dist_to_methano = np.linalg.norm(get_coords(conf, ester_carbons) - get_coords(conf, methano_idx))
        
        # Endo product: C=C bond is closer to the ethano-diester bridge than the methano bridge is.
        data['is_endo'] = dist_to_cc < dist_to_methano

    # --- Final Verdict ---
    verdict = ""
    for key in ['A', 'B', 'C', 'D']:
        data = options[key]
        if 'error' in data:
            if key == llm_answer:
                return f"Incorrect. The provided answer {key} could not be processed: {data['error']}"
            continue

        is_major_isomer = data['is_cis'] and data['is_anti'] and data['is_endo']
        
        if key == llm_answer:
            if is_major_isomer:
                verdict = "Correct"
            else:
                reasons = []
                if not data['is_cis']:
                    reasons.append("it has a trans-diester configuration, but the reaction requires a cis-diester")
                if not data['is_anti']:
                    reasons.append("it is a syn-adduct, but the sterically favored product is the anti-adduct")
                if not data['is_endo']:
                    reasons.append("it is an exo-adduct, but the kinetically favored product is the endo-adduct")
                return f"Incorrect. The provided answer {key} is wrong because {', and '.join(reasons)}."

    # This part will only be reached if the verdict on the LLM's answer was "Correct"
    # We double-check that no other option was also deemed correct.
    correct_isomers = [k for k, v in options.items() if 'error' not in v and v['is_cis'] and v['is_anti'] and v['is_endo']]
    if len(correct_isomers) > 1:
        return f"Check inconclusive. Multiple options ({', '.join(correct_isomers)}) satisfy all constraints. The question implies a single major product."
    
    return verdict

result = check_organic_synthesis_answer()
print(result)