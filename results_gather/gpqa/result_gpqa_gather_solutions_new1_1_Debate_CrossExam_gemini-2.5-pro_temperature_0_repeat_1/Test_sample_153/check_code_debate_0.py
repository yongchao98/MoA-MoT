def check_correctness():
    """
    This function checks the correctness of the proposed answer by verifying 
    each candidate molecule against the given spectral data constraints.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        return "Could not run the check. Please install the RDKit library using 'pip install rdkit'."

    # --- Define Candidates and Proposed Answer ---
    candidates = {
        'A': {'name': '2-chlorobenzoic acid', 'smiles': 'O=C(O)c1ccccc1Cl'},
        'B': {'name': 'Phenyl chloroformate', 'smiles': 'O=C(Cl)Oc1ccccc1'},
        'C': {'name': '3-Chloro-2-hydroxybenzaldehyde', 'smiles': 'O=Cc1c(O)cccc1Cl'},
        'D': {'name': '4-chlorobenzoic acid', 'smiles': 'O=C(O)c1ccc(Cl)cc1'}
    }
    proposed_answer_key = 'D'

    # --- Constraint Checking Functions ---

    def check_mass_spec(mol):
        """Checks for one Cl atom and a nominal mass of 156."""
        # 1. Check for exactly one chlorine atom
        has_one_cl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17) == 1
        if not has_one_cl:
            return False, "Does not contain exactly one chlorine atom."

        # 2. Check for nominal mass of 156 (with 35Cl)
        # Create a mutable copy to set the isotope for mass calculation
        rw_mol = Chem.RWMol(mol)
        cl_atom_found = False
        for atom in rw_mol.GetAtoms():
            if atom.GetAtomicNum() == 17:
                atom.SetIsotope(35)
                cl_atom_found = True
                break
        
        # Calculate mass using the 35Cl isotope
        mass = Descriptors.ExactMolWt(rw_mol)
        if not (155.9 < mass < 156.1):
            return False, f"Incorrect molecular weight. Calculated mass with 35Cl is {mass:.2f}, expected ~156."
            
        return True, "Passes MS check."

    def check_ir_spec(mol):
        """Checks for a carboxylic acid group based on IR data."""
        # Pattern for carboxylic acid: [C](=O)[O-H]
        cooh_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        if not mol.HasSubstructMatch(cooh_pattern):
            return False, "Lacks a carboxylic acid group, which is required by the broad 3500-2700 cm-1 and 1720 cm-1 IR peaks."
        return True, "Passes IR check."

    def check_nmr_spec(mol):
        """Checks for a para-disubstituted benzene ring."""
        # The 11.0 ppm signal is for the COOH proton, confirmed by the IR check.
        # The key is the aromatic pattern: two 2H doublets -> para-substitution.
        
        # Pattern for a 1,4-disubstituted benzene ring
        para_pattern = Chem.MolFromSmarts('c1ccc(c(c1)*)*') # * is any atom
        if mol.HasSubstructMatch(para_pattern):
            return True, "Passes NMR check (para-substituted)."
        
        # Check for other substitution patterns to provide a more specific reason for failure
        ortho_pattern = Chem.MolFromSmarts('c1cccc(c1*)*')
        if mol.HasSubstructMatch(ortho_pattern):
            return False, "Is ortho-substituted, which would give a complex NMR pattern, not two doublets."
            
        meta_pattern = Chem.MolFromSmarts('c1cc(cc(c1)*)*')
        if mol.HasSubstructMatch(meta_pattern):
            return False, "Is meta-substituted, which would give a complex NMR pattern, not two doublets."

        return False, "Does not have a para-disubstituted benzene ring required by the 2H+2H doublet pattern."

    # --- Run Checks and Evaluate ---
    passing_candidates = []
    results_log = {}

    for key, data in candidates.items():
        mol = Chem.MolFromSmiles(data['smiles'])
        if mol is None:
            results_log[key] = "Failed: Invalid SMILES string."
            continue

        ms_ok, ms_reason = check_mass_spec(mol)
        if not ms_ok:
            results_log[key] = f"Failed MS check: {ms_reason}"
            continue

        ir_ok, ir_reason = check_ir_spec(mol)
        if not ir_ok:
            results_log[key] = f"Failed IR check: {ir_reason}"
            continue
        
        nmr_ok, nmr_reason = check_nmr_spec(mol)
        if not nmr_ok:
            results_log[key] = f"Failed NMR check: {nmr_reason}"
            continue
        
        # If all checks pass
        results_log[key] = "Passed all checks."
        passing_candidates.append(key)

    # --- Final Verification ---
    if len(passing_candidates) == 0:
        return "Incorrect: No candidate satisfied all the spectral constraints."
    
    if len(passing_candidates) > 1:
        return f"Incorrect: The data is ambiguous as multiple candidates ({', '.join(passing_candidates)}) satisfied all constraints."

    # Exactly one candidate passed
    identified_key = passing_candidates[0]
    
    if identified_key == proposed_answer_key:
        return "Correct"
    else:
        failure_reason = results_log.get(proposed_answer_key, "the proposed answer was not evaluated due to an earlier failure.")
        return (f"Incorrect: The analysis shows that only candidate {identified_key} ({candidates[identified_key]['name']}) fits all the data. "
                f"The proposed answer {proposed_answer_key} was wrong. Reason: {failure_reason}")

# Execute the check and print the result
print(check_correctness())