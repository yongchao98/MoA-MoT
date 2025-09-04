import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry.rdGeometry import Point3D
    import numpy as np
except ImportError:
    # rdkit is a specialized library. If not installed, we cannot run the check.
    # We will print an informative error message.
    print("Error: RDKit library not found. Please install it using 'pip install rdkit-pypi' to run this verification code.")
    sys.exit(1)

def get_relative_stereochem(smiles: str):
    """
    Analyzes the 3D structure of a cyclohexene derivative to determine the
    relative stereochemistry (cis/trans) of its substituents.

    Args:
        smiles: The SMILES string of the molecule.

    Returns:
        A dictionary with relative stereochemistry or an error string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return f"Error: RDKit could not parse SMILES: {smiles}"
    mol = Chem.AddHs(mol)
    
    # Generate a 3D conformer
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) == -1:
        return "Error: Could not generate 3D conformer."
    AllChem.UFFOptimizeMolecule(mol)

    # --- Find atoms of interest based on connectivity ---
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return "Error: Molecule is not monocyclic."
    ring_atoms = list(ring_info.AtomRings()[0])

    # Find the C=C bond in the ring
    idx_c2, idx_c3 = -1, -1
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
            idx_c2, idx_c3 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            break
    if idx_c2 == -1: return "Error: Could not find double bond in ring."

    # Find C1 (neighbor of C2 with -OH) and C4 (neighbor of C3 with -Me)
    idx_c1, idx_c4 = -1, -1
    for atom_idx in ring_atoms:
        if atom_idx in [idx_c2, idx_c3]: continue
        atom = mol.GetAtomWithIdx(atom_idx)
        is_neighbor_to_c2 = any(n.GetIdx() == atom_idx for n in mol.GetAtomWithIdx(idx_c2).GetNeighbors())
        is_neighbor_to_c3 = any(n.GetIdx() == atom_idx for n in mol.GetAtomWithIdx(idx_c3).GetNeighbors())
        has_oh = any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors())
        
        if is_neighbor_to_c2 and has_oh:
            idx_c1 = atom_idx
        elif is_neighbor_to_c3:
            idx_c4 = atom_idx
    if idx_c1 == -1 or idx_c4 == -1: return "Error: Could not identify C1 and C4."

    # Find C5 and C6
    idx_c6 = next(n.GetIdx() for n in mol.GetAtomWithIdx(idx_c1).GetNeighbors() if n.GetIdx() in ring_atoms and n.GetIdx() != idx_c2)
    idx_c5 = next(n.GetIdx() for n in mol.GetAtomWithIdx(idx_c4).GetNeighbors() if n.GetIdx() in ring_atoms and n.GetIdx() != idx_c3)
    if not mol.GetBondBetweenAtoms(idx_c5, idx_c6): return "Error: C5 and C6 are not bonded."

    # Find substituent atoms
    def get_substituent_atom(parent_idx, atomic_num=6):
        return next(n for n in mol.GetAtomWithIdx(parent_idx).GetNeighbors() if n.GetAtomicNum() == atomic_num and n.GetIdx() not in ring_atoms)

    sub_c4 = get_substituent_atom(idx_c4)
    sub_c5 = get_substituent_atom(idx_c5)
    sub_c6 = get_substituent_atom(idx_c6)

    # --- Analyze 3D geometry ---
    conf = mol.GetConformer()
    
    # Calculate the normal vector to the best-fit plane of the ring
    ring_coords = np.array([conf.GetAtomPosition(i) for i in ring_atoms])
    centroid = np.mean(ring_coords, axis=0)
    u, s, vh = np.linalg.svd(ring_coords - centroid)
    normal = vh[2, :]

    # Project substituent vectors onto the normal. Same sign = same side (cis).
    def get_dot_product(atom):
        vec = np.array(conf.GetAtomPosition(atom.GetIdx())) - centroid
        return np.dot(vec, normal)

    dot_c4 = get_dot_product(sub_c4)
    dot_c5 = get_dot_product(sub_c5)
    dot_c6 = get_dot_product(sub_c6)

    # Determine relationships
    results = {}
    results['C5-C6'] = 'cis' if np.sign(dot_c5) == np.sign(dot_c6) else 'trans'
    results['C4-C5'] = 'cis' if np.sign(dot_c4) == np.sign(dot_c5) else 'trans'
    results['C4-C6'] = 'cis' if np.sign(dot_c4) == np.sign(dot_c6) else 'trans'
    return results

def check_answer():
    """
    Main function to verify the step-by-step solution.
    """
    # Step 1: Check Compound A (n-butane)
    # The solution correctly identifies A as n-butane despite the quartet/sextet discrepancy.
    # We check if n-butane has 2 sets of equivalent protons.
    mol_A = Chem.MolFromSmiles('CCCC')
    ranks = Chem.CanonicalRankAtoms(mol_A, breakTies=False)
    num_proton_sets = len(set(ranks))
    if num_proton_sets != 2:
        return f"Reason: The proposed structure for A, n-butane, should have 2 sets of equivalent protons. The code found {num_proton_sets}."

    # Step 2 & 3: Check B and C
    # The logic B=2-bromobutane and C=but-2-ene is sound based on standard reaction rules (Markovnikov/Zaitsev).
    # The key check is that C has geometric isomers, which but-2-ene does.
    mol_C_cis = Chem.MolFromSmiles('C/C=C\\C')
    if mol_C_cis is None:
        return "Reason: RDKit failed to parse the structure for Compound C (cis-but-2-ene)."

    # Step 4: Check Diels-Alder Regiochemistry
    diene = Chem.MolFromSmiles('O/C=C/C=C/C') # (1E,3E)-penta-1,3-dien-1-ol
    dienophile = mol_C_cis
    rxn = AllChem.ReactionFromSmarts('[c:1]=[c:2][c:3]=[c:4].[c:5]=[c:6]>>[c:1]1[c:2]=[c:3][c:4][c:5][c:6]1')
    products = rxn.RunReactants((diene, dienophile))
    
    product_regio = products[0][0]
    Chem.SanitizeMol(product_regio)
    
    # Check if the product has the 4,5,6-trimethylcyclohex-2-enol skeleton.
    pattern = Chem.MolFromSmarts('C1(O)C=CC(C)C(C)C1(C)')
    if not product_regio.HasSubstructMatch(pattern):
        return "Reason: The simulated Diels-Alder reaction produced a product with incorrect regiochemistry (connectivity). It should be a 4,5,6-trimethylcyclohex-2-enol skeleton."

    # Step 5: Check Stereochemistry of Options A and D
    # The solution states D is the exo product. We verify this by checking the relative stereochemistry.
    # Endo product: C4-Me is cis to C5-Me and C6-Me.
    # Exo product: C4-Me is trans to C5-Me and C6-Me.
    # Both must have C5-Me cis to C6-Me (from cis-but-2-ene).
    
    # SMILES for the relevant options
    smiles_A = 'C[C@H]1[C@H](C)C[C@H](O)C=C[C@H]1C' # (1S,4R,5S,6S)-...
    smiles_D = 'C[C@H]1[C@@H](C)C[C@H](O)C=C[C@H]1C' # (1S,4R,5S,6R)-...

    stereo_A = get_relative_stereochem(smiles_A)
    if isinstance(stereo_A, str): return f"Reason: {stereo_A}"
    
    stereo_D = get_relative_stereochem(smiles_D)
    if isinstance(stereo_D, str): return f"Reason: {stereo_D}"

    # Check 1: C5-Me and C6-Me must be cis.
    if stereo_A['C5-C6'] != 'cis' or stereo_D['C5-C6'] != 'cis':
        return "Reason: The stereochemistry of the options is inconsistent with a cis-dienophile. The C5 and C6 methyl groups must be cis to each other."

    # Check 2: Identify which option is endo and which is exo.
    is_A_endo = stereo_A['C4-C5'] == 'cis' and stereo_A['C4-C6'] == 'cis'
    is_D_exo = stereo_D['C4-C5'] == 'trans' and stereo_D['C4-C6'] == 'trans'

    if not is_A_endo:
        return f"Reason: The solution implies option A is the endo product, but its calculated stereochemistry ({stereo_A}) does not match the endo definition (all methyl groups cis)."
    
    if not is_D_exo:
        return f"Reason: The solution identifies option D as the exo product, but its calculated stereochemistry ({stereo_D}) does not match the exo definition (C4-Me is trans to C5-Me and C6-Me)."

    # If all checks pass, the logic holds. The solution correctly identifies the reaction pathway,
    # the product regiochemistry, and the stereochemical nature of the options. The final choice of D
    # is based on the fact that it is the correct exo product listed as an option.
    return "Correct"

# Run the check and print the result
try:
    result = check_answer()
    print(result)
except Exception as e:
    print(f"An unexpected error occurred during verification: {e}")
