from rdkit import Chem

def check_organic_synthesis_answer():
    """
    Checks the proposed final product of a multi-step synthesis against
    the chemical principles of each reaction step.
    """
    # The proposed final product is (2S,3R,4S,6S)-2-benzyl-4-hydroxy-6-methyl-3-phenylcyclohexan-1-one.
    # Its canonical isomeric SMILES from PubChem (CID 139031133) is used here.
    smiles = 'C[C@H]1CC[C@H]([C@H]([C@@H](C1=O)Cc2ccccc2)c3ccccc3)O'
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Incorrect: The chemical structure corresponding to the answer could not be parsed from its SMILES string."
        
        # This is crucial for RDKit to calculate and store R/S labels
        Chem.AssignAtomChiralTagsFromStructure(mol)

        # --- Define Patterns for Substructure Matching ---
        carbonyl_pattern = Chem.MolFromSmarts('[CX3](=O)')
        
        # --- Find Key Atoms and Positions ---
        carbonyl_match = mol.GetSubstructMatch(carbonyl_pattern)
        if not carbonyl_match:
            return "Incorrect: The final product must be a ketone, but no carbonyl group was found."
        
        carbonyl_c_idx = carbonyl_match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Find the two alpha-carbons adjacent to the carbonyl
        alpha_carbons = [n for n in carbonyl_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(alpha_carbons) != 2:
            return "Incorrect: The ketone is not part of a standard six-membered ring."

        # --- Constraint 1: Check Regiochemistry of Alkylations (Steps 2 & 3) ---
        # One alpha-carbon must have a benzyl group (C2), the other a methyl group (C6).
        c2_atom, c6_atom = None, None
        for alpha_c in alpha_carbons:
            has_benzyl = any("c1ccccc1C" in Chem.MolToSmarts(neighbor) for neighbor in alpha_c.GetNeighbors())
            has_methyl = any(neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1 for neighbor in alpha_c.GetNeighbors())
            
            if has_benzyl and not has_methyl:
                c2_atom = alpha_c
            elif has_methyl and not has_benzyl:
                c6_atom = alpha_c
            elif has_benzyl and has_methyl:
                 return "Incorrect: Methylation occurred at C2. LDA should deprotonate the less hindered C6 position under kinetic control."

        if c2_atom is None:
            return "Incorrect: The benzyl group from Step 2 is missing or at the wrong position. It should be at C2."
        if c6_atom is None:
            return "Incorrect: The methyl group from Step 3 is missing or at the wrong position. It should be at C6."

        # --- Constraint 2: Check Regiochemistry of Conjugate Addition and Hydroxyl Position ---
        # Find C3 (adjacent to C2) and C4 (adjacent to C3)
        c3_atom = next((n for n in c2_atom.GetNeighbors() if n.IsInRing() and n.GetIdx() != carbonyl_c_idx), None)
        if c3_atom is None: return "Incorrect: Ring structure is malformed around C2."
        
        has_phenyl = any(n.GetIsAromatic() for n in c3_atom.GetNeighbors())
        if not has_phenyl:
            return "Incorrect: The phenyl group from Step 2 is missing or at the wrong position. It should be at C3."

        c4_atom = next((n for n in c3_atom.GetNeighbors() if n.IsInRing() and n.GetIdx() != c2_atom.GetIdx()), None)
        if c4_atom is None: return "Incorrect: Ring structure is malformed around C3."

        has_hydroxyl = any(n.GetAtomicNum() == 8 for n in c4_atom.GetNeighbors())
        if not has_hydroxyl:
            return "Incorrect: The hydroxyl group is missing or at the wrong position. It should be at C4."

        # --- Constraint 3: Check Stereochemistry ---
        expected_stereo = {
            c2_atom.GetIdx(): 'S',
            c3_atom.GetIdx(): 'R',
            c4_atom.GetIdx(): 'S',
            c6_atom.GetIdx(): 'S'
        }
        
        mismatches = []
        for idx, expected_config in expected_stereo.items():
            atom = mol.GetAtomWithIdx(idx)
            if not atom.HasProp('_CIPCode') or atom.GetProp('_CIPCode') != expected_config:
                actual_config = atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else 'None'
                mismatches.append(f"C{idx+1} expected {expected_config}, found {actual_config}")

        if mismatches:
            return f"Incorrect: The stereochemistry is wrong. Mismatches: {'; '.join(mismatches)}."

        # If all checks pass, the answer is correct.
        return "Correct"

    except Exception as e:
        return f"An error occurred during verification: {e}"

# Run the check
result = check_organic_synthesis_answer()
print(result)