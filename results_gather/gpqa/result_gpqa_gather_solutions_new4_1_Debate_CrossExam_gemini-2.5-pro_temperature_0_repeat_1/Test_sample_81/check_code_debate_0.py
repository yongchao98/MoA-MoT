import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the stereochemistry
    of the provided options based on chemical principles.
    """
    
    # Data from the question
    smiles_options = {
        "A": "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        "B": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O",
        "C": "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O",
        "D": "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O"
    }
    
    names = {
        "A": "dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        "B": "dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        "C": "dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate",
        "D": "dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate"
    }

    llm_answer = "C"

    # --- Constraint 1: Check for cis-diester configuration ---
    # Maleic anhydride is a cis-dienophile, so the product must have cis-esters.
    # We check this from the IUPAC names given for C10 and C11.
    # (R,S) or (S,R) is cis. (R,R) or (S,S) is trans.
    cis_candidates = []
    for label, name in names.items():
        if "(10R,11R)" in name or "(10S,11S)" in name:
            pass  # This is a trans-diester
        else:
            cis_candidates.append(label)

    if llm_answer not in cis_candidates:
        return f"Incorrect. The product must have cis-diester groups because it originates from maleic anhydride. The IUPAC name for option {llm_answer} indicates a trans configuration for the ester groups at C10 and C11, which is incorrect."

    # --- Constraint 2: Check for anti-bridge configuration ---
    # The major product should be the 'anti' isomer due to sterics.
    # We determine this geometrically from the 3D structure.
    
    def get_bridge_orientation(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, randomSeed=42) == -1: return "conformer_failed"
        if AllChem.UFFOptimizeMolecule(mol) == -1: return "optimization_failed"
        
        conf = mol.GetConformer()
        rings = mol.GetRingInfo().AtomRings()
        
        # Find central 4-membered ring and its centroid
        central_ring_atoms = [r for r in rings if len(r) == 4][0]
        centroid = np.mean([conf.GetAtomPosition(i) for i in central_ring_atoms], axis=0)
        
        # Find methano bridge atom (CH2 in a 5-ring, not in central ring)
        methano_bridge_atom_idx = -1
        for r in rings:
            if len(r) == 5:
                unique_atoms = set(r) - set(central_ring_atoms)
                if len(unique_atoms) == 1:
                    methano_bridge_atom_idx = list(unique_atoms)[0]
                    break
        
        # Find ethano bridge atoms (2x CH in a 6-ring, not in central ring)
        ethano_bridge_atoms = []
        for r in rings:
            if len(r) == 6:
                unique_atoms = set(r) - set(central_ring_atoms)
                # The ethano bridge is part of the bicyclo[2.2.2]octene system
                if len(unique_atoms) == 2:
                    ethano_bridge_atoms = list(unique_atoms)
                    break

        if methano_bridge_atom_idx == -1 or len(ethano_bridge_atoms) != 2:
            return "bridge_id_failed"

        # Create vectors from centroid to bridges
        vec_methano = conf.GetAtomPosition(methano_bridge_atom_idx) - centroid
        ethano_midpoint = (conf.GetAtomPosition(ethano_bridge_atoms[0]) + conf.GetAtomPosition(ethano_bridge_atoms[1])) / 2.0
        vec_ethano = ethano_midpoint - centroid
        
        # Calculate angle between vectors
        dot_product = np.dot(vec_methano, vec_ethano) / (np.linalg.norm(vec_methano) * np.linalg.norm(vec_ethano))
        angle_deg = np.degrees(np.arccos(np.clip(dot_product, -1.0, 1.0)))
        
        return "anti" if angle_deg > 120 else "syn"

    orientations = {label: get_bridge_orientation(smiles_options[label]) for label in cis_candidates}
    
    anti_isomers = [label for label, orient in orientations.items() if orient == "anti"]
    
    if len(anti_isomers) != 1:
        return f"Error: Could not uniquely identify the anti-isomer among the valid candidates. Orientations found: {orientations}"
        
    predicted_major_product = anti_isomers[0]

    if llm_answer == predicted_major_product:
        return "Correct"
    else:
        return f"Incorrect. The major product should have an 'anti' configuration of the two main bridges due to sterically favored addition. The code identifies structure {predicted_major_product} as the anti-isomer, but the provided answer was {llm_answer}. The calculated bridge orientations for the cis-isomers were: {orientations}."

# Execute the check
result = check_correctness()
print(result)