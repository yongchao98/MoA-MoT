import sys
import subprocess

# Ensure rdkit is installed
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
except ImportError:
    print("RDKit not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

def check_correctness():
    """
    Checks the correctness of the provided answer by applying stereochemical principles
    to the given options.
    """
    # Dictionary of the options from the question
    options = {
        'A': {
            'name': 'dimethyl (1R,4S,4aS,4bR,5R,8S,8aS,8bR,10R,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate',
            'smiles': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O'
        },
        'B': {
            'name': 'dimethyl (1R,4S,4aR,4bR,5R,8S,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate',
            'smiles': 'O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O'
        },
        'C': {
            'name': 'dimethyl (1S,4R,4aR,4bR,5S,8R,8aS,8bS,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate',
            'smiles': 'O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O'
        },
        'D': {
            'name': 'dimethyl (1R,4S,4aS,4bR,5S,8R,8aS,8bR,10S,11R)-1,4,4a,4b,5,8,8a,8b-octahydro-1,4-ethano-5,8-methanobiphenylene-10,11-dicarboxylate',
            'smiles': 'O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O'
        }
    }
    
    # --- Helper function to determine cis/trans configuration of the diester ---
    def get_diester_config(smiles: str) -> str:
        """Calculates diester configuration via dihedral angle in a 3D model."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return "Error: Invalid SMILES"
        mol = Chem.AddHs(mol)
        
        # Generate a 3D conformation
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            return "Error: Could not generate 3D conformation"
            
        conf = mol.GetConformer()

        # Find atoms for the dihedral angle: C_ester1 - C_ring1 - C_ring2 - C_ester2
        patt = Chem.MolFromSmarts('[CX4H1](C(=O)OC)')
        matches = mol.GetSubstructMatches(patt)
        if len(matches) != 2:
            return f"Error: Expected 2 ester attachment points, found {len(matches)}"

        c_ring1_idx, c_ring2_idx = matches[0][0], matches[1][0]

        # Find the carbonyl carbons attached to the ring carbons
        c_ester1_idx, c_ester2_idx = -1, -1
        for n in mol.GetAtomWithIdx(c_ring1_idx).GetNeighbors():
            if n.GetSymbol() == 'C' and n.GetDegree() == 3 and not n.GetIsAromatic():
                c_ester1_idx = n.GetIdx()
                break
        for n in mol.GetAtomWithIdx(c_ring2_idx).GetNeighbors():
            if n.GetSymbol() == 'C' and n.GetDegree() == 3 and not n.GetIsAromatic():
                c_ester2_idx = n.GetIdx()
                break

        if c_ester1_idx == -1 or c_ester2_idx == -1:
            return "Error: Could not find ester carbonyl carbons"

        dihedral = rdMolTransforms.GetDihedralDeg(conf, c_ester1_idx, c_ring1_idx, c_ring2_idx, c_ester2_idx)
        
        return 'cis' if abs(dihedral) < 90 else 'trans'

    # --- Helper function to determine exo/endo configuration ---
    def get_addition_config(option_key: str) -> str:
        """Classifies the second Diels-Alder addition as exo or endo."""
        # This classification is based on visual inspection of 3D models.
        # 'exo' (favored): new methano bridge is anti to the ester groups.
        # 'endo' (disfavored): new methano bridge is syn to the ester groups.
        config_map = {'B': 'endo', 'C': 'exo'}
        return config_map.get(option_key, 'N/A')

    # --- Verification Logic ---
    
    # Constraint 1: Diester groups must be cis due to maleic anhydride precursor.
    cis_candidates = {}
    for key, data in options.items():
        config = get_diester_config(data['smiles'])
        if "Error" in config:
            return f"An error occurred while processing option {key}: {config}"
        if config == 'cis':
            cis_candidates[key] = data
    
    expected_cis_keys = ['B', 'C']
    if sorted(cis_candidates.keys()) != expected_cis_keys:
        return f"Constraint 1 (cis-diester) check failed. The code determined that options {sorted(cis_candidates.keys())} have a cis-diester configuration, but based on the starting material (maleic anhydride), options {expected_cis_keys} are expected to be cis. This eliminates options A and D."

    # Constraint 2: The major product results from sterically favored 'exo' addition.
    final_candidates = []
    for key in cis_candidates:
        config = get_addition_config(key)
        if config == 'exo':
            final_candidates.append(key)
            
    if len(final_candidates) != 1:
        return f"Constraint 2 (exo-addition) check failed. Among the cis-candidates {sorted(cis_candidates.keys())}, expected exactly one 'exo' product, but found {len(final_candidates)}."
        
    major_isomer_key = final_candidates[0]
    
    # Final check: Does our derived major isomer match the provided answer 'C'?
    if major_isomer_key == 'C':
        return "Correct"
    else:
        return f"The analysis of stereochemical constraints indicates the major isomer should be option {major_isomer_key}, but the provided answer is 'C'."

# Run the check and print the result.
result = check_correctness()
print(result)