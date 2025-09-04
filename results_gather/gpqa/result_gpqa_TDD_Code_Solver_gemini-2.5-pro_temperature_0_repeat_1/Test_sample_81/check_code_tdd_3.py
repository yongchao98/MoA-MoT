import sys

try:
    from rdkit import Chem
except ImportError:
    # This check is for the user's environment.
    # If rdkit is not installed, the check cannot be performed.
    print("Constraint Not Satisfied: RDKit library not found. Please install it to run this check, e.g., 'pip install rdkit-pypi'")
    sys.exit(1)

def check_correctness():
    '''
    This function checks the correctness of the provided answer by verifying the 
    stereochemical claims made in the explanation against the provided SMILES strings.
    '''
    # Data from the question and options
    data = {
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

    # The provided answer's logical claims to be verified
    explanation_claims = {
        "cis_diester": {
            "trans": ["A"],
            "cis": ["B", "C", "D"]
        }
    }

    def find_ester_bearing_carbons(mol):
        '''Finds the two chiral carbons bonded to the ester groups.'''
        ester_bearing_carbons = []
        # SMARTS pattern for a methyl ester group's carbonyl carbon
        patt = Chem.MolFromSmarts('[CX3](=O)[O][CX4H3]')
        matches = mol.GetSubstructMatches(patt)
        
        if len(matches) != 2:
            return None
        
        carbonyl_indices = [m[0] for m in matches]
        
        for c_idx in carbonyl_indices:
            carbonyl_atom = mol.GetAtomWithIdx(c_idx)
            for neighbor in carbonyl_atom.GetNeighbors():
                # The neighbor that is a chiral carbon is the one we want
                if neighbor.GetAtomicNum() == 6 and neighbor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    ester_bearing_carbons.append(neighbor)
        
        # Use a set to count unique atoms found
        if len(set(atom.GetIdx() for atom in ester_bearing_carbons)) != 2:
            return None
            
        return ester_bearing_carbons

    def get_diester_stereochem(smiles):
        '''
        Determines if the diester is cis or trans based on the R/S configuration
        of the carbons they are attached to.
        '''
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return "invalid_smiles"
        
        carbons = find_ester_bearing_carbons(mol)
        if not carbons:
            return "ester_carbons_not_found"
        
        atom1, atom2 = carbons
        
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        
        try:
            cip1 = atom1.GetProp('_CIPCode')
            cip2 = atom2.GetProp('_CIPCode')
        except KeyError:
            return "cip_assignment_failed"

        # For a 1,2-disubstituted system on this 5-membered ring,
        # different CIP codes (R,S) mean cis, same codes (R,R or S,S) mean trans.
        return 'cis' if cip1 != cip2 else 'trans'

    # --- Verification Steps ---

    # 1. Verify the cis/trans claim from the explanation
    # Claim: Maleic anhydride leads to a cis-diester, eliminating trans-isomers.
    for option, details in data.items():
        stereochem = get_diester_stereochem(details['smiles'])
        
        is_trans_claim = option in explanation_claims['cis_diester']['trans']
        is_cis_claim = option in explanation_claims['cis_diester']['cis']

        if is_trans_claim and stereochem != 'trans':
            return f"Incorrect: The explanation states that option {option} should be a trans-diester, but the provided SMILES corresponds to a '{stereochem}' configuration. The reasoning is based on a false premise."
        if is_cis_claim and stereochem != 'cis':
            return f"Incorrect: The explanation states that option {option} should be a cis-diester, but the provided SMILES corresponds to a '{stereochem}' configuration. The reasoning is based on a false premise."

    # 2. Verify that B, C, and D are distinct stereoisomers as claimed
    # Claim: B, C, and D are different isomers (endo, exo from syn-precursor, exo from anti-precursor).
    mols = {opt: Chem.MolFromSmiles(data[opt]['smiles']) for opt in ['B', 'C', 'D']}
    for opt, mol in mols.items():
        if not mol:
            return f"Incorrect: The SMILES string for option {opt} is invalid and cannot be parsed."

    # Canonical SMILES is a unique representation for a molecule, allowing for comparison.
    canon_smiles = {opt: Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) for opt, mol in mols.items()}

    if canon_smiles['B'] == canon_smiles['D']:
        return "Incorrect: The explanation claims B (endo) and D (exo) are different isomers, but their SMILES represent the same molecule."
    
    if canon_smiles['C'] == canon_smiles['D']:
        return "Incorrect: The explanation claims C (from syn-precursor) and D (from anti-precursor) are different isomers, but their SMILES represent the same molecule."

    if canon_smiles['B'] == canon_smiles['C']:
        return "Incorrect: The explanation presents options B and C as distinct isomers to be eliminated separately, but their SMILES represent the same molecule."

    # 3. Final Conclusion
    # The code has verified the structural properties asserted in the explanation:
    # - A is trans, while B, C, D are cis. This correctly eliminates A.
    # - B, C, and D are all distinct stereoisomers. This means the subsequent elimination steps are based on real structural differences.
    # The chemical reasoning itself (e.g., "exo-addition is favored") is a principle of organic chemistry that this code does not verify.
    # However, since the structural premises of the reasoning are correct, the logical conclusion is valid.

    return "Correct"

# Run the check
try:
    result = check_correctness()
    print(result)
except Exception as e:
    print(f"An unexpected error occurred during the check: {e}")