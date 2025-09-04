def check_correctness_of_chemistry_answer():
    """
    This function checks the correctness of the LLM's answer for a chemical identification problem.
    It uses the RDKit library to analyze the properties of the proposed molecules against the given spectral data.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError:
        return "Could not run the check because the 'rdkit' library is not installed. Please install it using 'pip install rdkit-pypi'."

    # --- Problem Data ---
    # The answer provided by the LLM
    llm_answer_key = "B"

    # Candidate molecules and their SMILES representations
    candidates = {
        "A": {"name": "Phenyl chloroformate", "smiles": "O=C(O[c]1[cH][cH][cH][cH][cH]1)Cl"},
        "B": {"name": "4-chlorobenzoic acid", "smiles": "O=C(O)c1ccc(Cl)cc1"},
        "C": {"name": "2-chlorobenzoic acid", "smiles": "O=C(O)c1ccccc1Cl"},
        "D": {"name": "3-Chloro-2-hydroxybenzaldehyde", "smiles": "O=Cc1c(O)ccc(Cl)c1"}
    }

    # Constraints derived from the spectral data
    constraints = {
        "m_plus_mass": 156.0,
        "m_plus_2_mass": 158.0,
        "has_chlorine": True,
        "has_carboxylic_acid": True,
        "num_aromatic_protons": 4,
        "aromatic_proton_symmetry": [2, 2] # Expects two groups of two equivalent protons (para-substitution)
    }

    # --- Verification Logic ---
    results = {}
    for key, data in candidates.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if not mol:
            results[key] = f"Invalid SMILES string for {data['name']}."
            continue
        
        mol = Chem.AddHs(mol)
        
        # 1. Mass Spec Check
        cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
        if constraints["has_chlorine"] and len(cl_atoms) != 1:
            results[key] = f"Failed MS check: Expected 1 chlorine atom, but found {len(cl_atoms)}."
            continue
        
        # Calculate M+ peak (with 35Cl, the most abundant isotope)
        m_plus_mass = rdMolDescriptors.CalcExactMolWt(mol)
        if not (constraints["m_plus_mass"] - 0.1 < m_plus_mass < constraints["m_plus_mass"] + 0.1):
            results[key] = f"Failed MS check: Calculated M+ mass is {m_plus_mass:.2f}, expected ~{constraints['m_plus_mass']}."
            continue

        # 2. IR Check (Functional Groups)
        carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
        has_cooh = mol.HasSubstructMatch(carboxylic_acid_pattern)
        if constraints["has_carboxylic_acid"] and not has_cooh:
            results[key] = "Failed IR check: Does not contain a carboxylic acid group, which is indicated by the broad 3500-2700 and 1720 cm-1 peaks."
            continue

        # 3. 1H NMR Check
        aromatic_protons = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic() and atom.GetAtomicNum() == 1]
        if len(aromatic_protons) != constraints["num_aromatic_protons"]:
            results[key] = f"Failed NMR check: Expected {constraints['num_aromatic_protons']} aromatic protons, but found {len(aromatic_protons)}."
            continue
            
        # Check symmetry of aromatic protons
        symmetry_classes = Chem.rdmolops.GetSymmetryClasses(mol)
        aromatic_proton_indices = [p.GetIdx() for p in aromatic_protons]
        
        aromatic_proton_symmetry_groups = {}
        for idx in aromatic_proton_indices:
            class_id = symmetry_classes[idx]
            aromatic_proton_symmetry_groups[class_id] = aromatic_proton_symmetry_groups.get(class_id, 0) + 1
        
        signal_integrations = sorted(list(aromatic_proton_symmetry_groups.values()))
        
        if signal_integrations != constraints["aromatic_proton_symmetry"]:
            results[key] = f"Failed NMR check: The aromatic proton symmetry is {signal_integrations}, which does not match the expected {constraints['aromatic_proton_symmetry']} pattern for two doublets of 2H each (para-substitution)."
            continue

        # If all checks pass for a candidate
        results[key] = "Passed all checks."

    # --- Final Verdict ---
    llm_result = results.get(llm_answer_key)
    
    if llm_result == "Passed all checks.":
        # Verify that no other candidate also passed
        other_passed = [k for k, v in results.items() if v == "Passed all checks." and k != llm_answer_key]
        if not other_passed:
            return "Correct"
        else:
            return f"The LLM's answer {llm_answer_key} is plausible, but the check also found candidate(s) {other_passed} passed all checks. The constraints may need refinement."
    else:
        return f"The LLM's answer {llm_answer_key} ({candidates[llm_answer_key]['name']}) is incorrect. Reason: {llm_result}"

# Execute the check and print the result
print(check_correctness_of_chemistry_answer())