def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's answer regarding optical isomerism.
    It uses the RDKit library to analyze the molecular structures.
    """
    # Try to import rdkit. If it's not installed, return an informative error.
    try:
        from rdkit import Chem
        from rdkit.Chem.rdmolops import FindMolChiralCenters
    except ImportError:
        return ("Could not run check: The 'rdkit' library is required. "
                "Please install it, for example, using 'pip install rdkit-pypi'.")

    def check_biphenyl_atropisomerism(mol):
        """
        A specific check for atropisomerism in biphenyls.
        Returns True if all four ortho positions are substituted by non-hydrogen atoms,
        which is a strong indicator of hindered rotation leading to chirality.
        """
        # Pattern for a single bond connecting two aromatic atoms
        patt = Chem.MolFromSmarts('[a;r]-[a;r]')
        matches = mol.GetSubstructMatches(patt)
        
        for match in matches:
            atom1_idx, atom2_idx = match
            atom1 = mol.GetAtomWithIdx(atom1_idx)
            atom2 = mol.GetAtomWithIdx(atom2_idx)
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)

            # Ensure it's a single bond between two different aromatic rings
            # (not part of a fused ring system like naphthalene)
            if not (bond.GetBondType() == Chem.BondType.SINGLE and \
                    atom1.GetRingInfo().NumAtomRings(atom1_idx) == 1 and \
                    atom2.GetRingInfo().NumAtomRings(atom2_idx) == 1):
                continue

            # Find ortho positions for each atom in the bond
            ortho_atoms = []
            for atom_in_bond, other_atom_idx in [(atom1, atom2_idx), (atom2, atom1_idx)]:
                for neighbor in atom_in_bond.GetNeighbors():
                    # An ortho neighbor is in the same ring but is not the other atom from the central bond
                    if neighbor.GetIdx() != other_atom_idx and neighbor.GetIsAromatic():
                        ortho_atoms.append(neighbor)
            
            # A valid biphenyl link should have 4 ortho atoms (2 for each ring)
            if len(ortho_atoms) != 4:
                continue

            # Check if all four ortho positions are substituted with non-hydrogen atoms.
            # An aromatic carbon's degree is 2 if unsubstituted, >2 if substituted.
            # We use GetTotalDegree() on the H-added molecule for accuracy.
            all_ortho_substituted = all(atom.GetTotalDegree() > 2 for atom in ortho_atoms)

            if all_ortho_substituted:
                return True # Found a potential atropisomeric axis
        return False

    def check_optical_isomerism(smiles_string, name):
        """
        Checks a molecule for common sources of optical isomerism.
        Returns a tuple (is_chiral, reason).
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return (False, f"Could not parse SMILES for {name}")

        # Add explicit hydrogens to get correct degree counts and find chiral centers
        mol_with_hs = Chem.AddHs(mol)

        # 1. Check for chiral centers (stereocenters)
        # includeUnassigned=True finds potential chiral centers even if not specified R/S
        chiral_centers = FindMolChiralCenters(mol_with_hs, includeUnassigned=True)
        if chiral_centers:
            return (True, f"has {len(chiral_centers)} chiral center(s)")

        # 2. Special check for biphenyl atropisomerism
        if "biphenyl" in name:
            if check_biphenyl_atropisomerism(mol_with_hs):
                return (True, "is a tetra-ortho-substituted biphenyl, which exhibits atropisomerism")

        # If no chiral centers or other known sources are found, assume achiral.
        return (False, "lacks chiral centers and other common sources of chirality like atropisomerism")

    # --- Main execution ---
    # Define the compounds from the question
    compounds = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "smiles": "COC(=O)c1c(cccc1[N+](=O)[O-])-c2c(C(=O)OC)cccc2[N+](=O)[O-]"
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "smiles": "CC(O)C(=O)OC"
        },
        3: {
            "name": "benzophenone",
            "smiles": "O=C(c1ccccc1)c2ccccc2"
        },
        4: {
            "name": "dimethyl fumarate",
            "smiles": "COC(=O)/C=C/C(=O)OC"
        }
    }

    # The LLM's answer identifies compounds 1 and 2 as optically active, corresponding to option C.
    llm_chiral_indices = {1, 2}
    llm_option = "C"

    # Perform the checks using the code
    code_chiral_indices = set()
    analysis_details = {}
    for key, data in compounds.items():
        is_chiral, reason = check_optical_isomerism(data["smiles"], data["name"])
        if is_chiral:
            code_chiral_indices.add(key)
        analysis_details[key] = f"Compound {key} ({data['name']}) {reason}."

    # Compare the code's findings with the LLM's answer
    if llm_chiral_indices == code_chiral_indices:
        # The set of compounds is correct. Now check if the option letter is correct.
        # Option C is "1 and 2".
        if llm_option == "C":
            return "Correct"
        else:
            return (f"Incorrect. The reasoning is correct, but the final option letter is wrong. "
                    f"The LLM correctly identified that compounds {sorted(list(code_chiral_indices))} show optical isomerism, "
                    f"but chose option {llm_option} instead of C.")
    else:
        reasons = []
        wrongly_included = sorted(list(llm_chiral_indices - code_chiral_indices))
        if wrongly_included:
            reasons.append(f"The LLM incorrectly claims compound(s) {wrongly_included} show optical isomerism. "
                           f"Analysis: {analysis_details[wrongly_included[0]]}")
        
        missed = sorted(list(code_chiral_indices - llm_chiral_indices))
        if missed:
            reasons.append(f"The LLM missed that compound(s) {missed} show optical isomerism. "
                           f"Analysis: {analysis_details[missed[0]]}")
        
        return f"Incorrect. The LLM's answer is wrong. {' '.join(reasons)}"

# Execute the check and print the result.
result = check_correctness_of_llm_answer()
print(result)