import sys
import io
from contextlib import redirect_stdout

# Suppress RDKit warnings if possible, as they are not relevant to the final output.
try:
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.warning')
except ImportError:
    pass

def check_correctness():
    """
    Checks the correctness of the given answer about a multi-step organic synthesis.
    The check relies on the RDKit library for cheminformatics analysis.
    It verifies each step of the reasoning provided in the answer by applying
    stereochemical and geometric constraints to the candidate molecules.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolTransforms
    except ImportError:
        return "Skipping check: RDKit library is not installed. This check requires RDKit."

    # Helper function for geometric calculations
    def get_centroid(mol, indices):
        conf = mol.GetConformer()
        if not indices:
            return None
        centroid = Chem.rdGeometry.Point3D(0, 0, 0)
        for idx in indices:
            centroid += conf.GetAtomPosition(idx)
        return centroid / len(indices)

    def get_vector_angle(v1, v2):
        # RDKit's AngleTo returns radians, convert to degrees
        return v1.AngleTo(v2) * 180.0 / 3.14159

    # Check 1: Cis/Trans configuration of the diester groups
    def check_cis_trans_esters(smiles):
        """
        Determines if the two ester groups are cis or trans using the dihedral angle
        between the ester substituents. A small dihedral (<90) indicates a cis relationship.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return "Error: Invalid SMILES"
        
        ester_pattern = Chem.MolFromSmarts('[C&H1]([C](=O)O[CH3])')
        matches = mol.GetSubstructMatches(ester_pattern)
        if len(matches) != 2:
            return "Error: Could not find two ester attachment points."
        
        c1_idx, c2_idx = matches[0][0], matches[1][0]
        
        mol = Chem.AddHs(mol)
        # Use a reproducible random seed for consistent conformer generation
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        
        # Find the carbonyl carbon of each ester group to define the dihedral
        c1_ester_C_idx = -1
        for n in c1.GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetTotalNumHs() == 0: # Carbonyl carbon
                c1_ester_C_idx = n.GetIdx()
                break
                
        c2_ester_C_idx = -1
        for n in c2.GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetTotalNumHs() == 0: # Carbonyl carbon
                c2_ester_C_idx = n.GetIdx()
                break
                
        if c1_ester_C_idx == -1 or c2_ester_C_idx == -1:
            return "Error: Could not find ester carbonyl carbons."

        try:
            # Dihedral angle: (ester_C1) - c1 - c2 - (ester_C2)
            dihedral = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(), c1_ester_C_idx, c1_idx, c2_idx, c2_ester_C_idx)
        except Exception:
            return "Error: Could not calculate dihedral angle."
        
        return "cis" if abs(dihedral) < 90 else "trans"

    # Check 2: Facial Selectivity (syn/anti)
    def check_syn_anti(mol):
        """
        Checks if the adduct is syn or anti.
        anti: methano bridge and ethano-diester bridge are on opposite sides of the molecule.
        syn: methano bridge and ethano-diester bridge are on the same side.
        This is determined by the angle between vectors from the molecule's center to each bridge.
        """
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        conf = mol_h.GetConformer()

        ester_pattern = Chem.MolFromSmarts('[C&H1]([C](=O)O[CH3])')
        ester_matches = mol_h.GetSubstructMatches(ester_pattern)
        ester_attachment_indices = [m[0] for m in ester_matches]
        if len(ester_attachment_indices) != 2: return "Error: ester bridge not found"
        ester_bridge_centroid = get_centroid(mol_h, ester_attachment_indices)

        methano_bridge_idx = -1
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 2:
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 2:
                    n1, n2 = neighbors
                    if n1.GetDegree() >= 3 and n2.GetDegree() >= 3 and mol_h.GetRingInfo().IsAtomInRingOfSize(atom.GetIdx(), 5):
                        methano_bridge_idx = atom.GetIdx()
                        break
        if methano_bridge_idx == -1: return "Error: methano bridge not found"
        methano_coords = conf.GetAtomPosition(methano_bridge_idx)

        all_indices = list(range(mol_h.GetNumAtoms()))
        mol_centroid = get_centroid(mol_h, all_indices)

        vec_ester = ester_bridge_centroid - mol_centroid
        vec_methano = methano_coords - mol_centroid
        
        angle = get_vector_angle(vec_ester, vec_methano)

        if angle > 120:
            return "anti"
        elif angle < 60:
            return "syn"
        else:
            return f"intermediate ({angle:.1f} deg)"

    # Check 3: Orientation (endo/exo)
    def check_endo_exo(mol):
        """
        Checks if the adduct is endo or exo.
        endo: The bulky substituent (diester bridge) is on the same side as the norbornene C=C bond.
        exo: The bulky substituent is on the same side as the norbornene methano (-CH2-) bridge.
        """
        mol_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_h, randomSeed=42)
        conf = mol_h.GetConformer()
        ri = mol_h.GetRingInfo()

        ester_pattern = Chem.MolFromSmarts('[C&H1]([C](=O)O[CH3])')
        ester_matches = mol_h.GetSubstructMatches(ester_pattern)
        ester_attachment_indices = [m[0] for m in ester_matches]
        if len(ester_attachment_indices) != 2: return "Error: ester bridge not found"
        ester_bridge_centroid = get_centroid(mol_h, ester_attachment_indices)

        methano_bridge_idx = -1
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 2:
                if ri.IsAtomInRingOfSize(atom.GetIdx(), 5) and not ri.IsAtomInRingOfSize(atom.GetIdx(), 6):
                     methano_bridge_idx = atom.GetIdx()
                     break
        if methano_bridge_idx == -1: return "Error: methano bridge not found"
        methano_coords = conf.GetAtomPosition(methano_bridge_idx)

        db_pattern = Chem.MolFromSmarts('C=C')
        db_match = mol_h.GetSubstructMatch(db_pattern)
        if not db_match: return "Error: double bond not found"
        db_centroid = get_centroid(mol_h, db_match)

        bridgehead_indices = []
        methano_atom = mol_h.GetAtomWithIdx(methano_bridge_idx)
        n1, n2 = methano_atom.GetNeighbors()
        for atom in list(n1.GetNeighbors()) + list(n2.GetNeighbors()):
            idx = atom.GetIdx()
            if idx not in [n1.GetIdx(), n2.GetIdx(), methano_bridge_idx] and ri.NumAtomRings(idx) >= 3:
                if idx not in bridgehead_indices:
                    bridgehead_indices.append(idx)
        if len(bridgehead_indices) != 2: return "Error: bridgeheads not found"
        bridgehead_centroid = get_centroid(mol_h, bridgehead_indices)

        vec_ester = ester_bridge_centroid - bridgehead_centroid
        vec_methano = methano_coords - bridgehead_centroid
        vec_db = db_centroid - bridgehead_centroid

        angle_ester_methano = get_vector_angle(vec_ester, vec_methano)
        angle_ester_db = get_vector_angle(vec_ester, vec_db)

        return "endo" if angle_ester_db < angle_ester_methano else "exo"

    # --- Main Checking Logic ---
    smiles_dict = {
        'A': "O=C(OC)[C@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@@H]5C=C[C@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'B': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@@H]4[C@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@H]1C(OC)=O",
        'C': "O=C(OC)[C@@H]1[C@@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@H]5C=C[C@@H]4C5)[C@@H]2[C@@H]1C(OC)=O",
        'D': "O=C(OC)[C@@H]1[C@H](C=C2)[C@@H]3[C@@H]([C@H]4[C@@H]3[C@@H]5C=C[C@H]4C5)[C@H]2[C@@H]1C(OC)=O"
    }
    
    # Step 1: Check cis/trans configuration.
    # Reaction with maleic anhydride (a cis-dienophile) requires the final diester to be cis.
    cis_trans_results = {name: check_cis_trans_esters(s) for name, s in smiles_dict.items()}
    
    if cis_trans_results.get('B') != 'trans':
        return f"Incorrect. The reasoning states candidate B is trans, but the code found it is {cis_trans_results.get('B')}. The first elimination step is based on this premise, which appears to be flawed if the code is correct."
    
    for name in ['A', 'C', 'D']:
        if cis_trans_results.get(name) != 'cis':
            return f"Incorrect. Candidate {name} should have a cis-diester configuration but was found to be {cis_trans_results.get(name)}."

    # Step 2: Check facial selectivity (syn/anti).
    # The major product should be the 'anti' adduct due to steric hindrance.
    mols = {name: Chem.MolFromSmiles(s) for name, s in smiles_dict.items()}
    syn_anti_results = {name: check_syn_anti(mols[name]) for name in ['A', 'C', 'D']}

    if syn_anti_results.get('A') != 'anti':
        return f"Incorrect. The reasoning states candidate A is the 'anti' adduct, but the code identified it as {syn_anti_results.get('A')}."
    if syn_anti_results.get('C') != 'syn':
        return f"Incorrect. The reasoning states candidate C is a 'syn' adduct, but the code identified it as {syn_anti_results.get('C')}."
    if syn_anti_results.get('D') != 'syn':
        return f"Incorrect. The reasoning states candidate D is a 'syn' adduct, but the code identified it as {syn_anti_results.get('D')}."

    # Step 3: Check orientation (endo/exo).
    # The major product should be the kinetic 'endo' product.
    endo_exo_result = check_endo_exo(mols['A'])

    if endo_exo_result != 'endo':
        return f"Incorrect. The reasoning states the final product A should be the 'endo' adduct, but the code identified it as {endo_exo_result}."

    # If all checks pass, the reasoning is sound and the answer is correct.
    return "Correct"

# Run the check and print the result.
# The output will be "Correct" if the answer and reasoning are validated,
# or an error message explaining the discrepancy.
result = check_correctness()
print(result)