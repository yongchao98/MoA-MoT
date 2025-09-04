def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the chirality of each compound.
    This function requires the RDKit library.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Could not perform check: RDKit library is not installed. Please run 'pip install rdkit'."

    # Define the compounds with their names and SMILES representations
    compounds = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "smiles": "COC(=O)c1c(cccc1[N+](=O)[O-])c2c(C(=O)OC)cccc2[N+](=O)[O-]"
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "smiles": "CC(C(=O)OC)O"
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

    # The LLM's answer is B, which corresponds to compounds 1 and 2 being optically active.
    llm_active_indices = [1, 2]

    # Our analysis results will be stored here
    code_active_indices = []
    reasons = {}

    # Analyze each compound
    for idx, data in compounds.items():
        is_active = False
        reason = ""
        
        # Special handling for atropisomerism (axial chirality)
        if idx == 1:
            # This biphenyl has bulky ortho-substituents (-NO2 and -COOCH3)
            # that hinder rotation around the central C-C bond. This creates
            # stable, non-superimposable mirror images (enantiomers).
            # This is a classic case of atropisomerism.
            is_active = True
            reason = "Shows optical isomerism due to atropisomerism (axial chirality)."
        else:
            mol = Chem.MolFromSmiles(data["smiles"])
            if not mol:
                reasons[idx] = f"Failed to parse SMILES for {data['name']}"
                continue

            # Find potential tetrahedral chiral centers
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            if len(chiral_centers) > 0:
                # Compound 2 has a chiral carbon
                is_active = True
                reason = f"Is optically active because it has a chiral center at atom {chiral_centers[0][0]}."
            else:
                # Compounds 3 and 4 are achiral
                is_active = False
                if idx == 3:
                    reason = "Is achiral because it is a symmetrical molecule with a plane of symmetry."
                elif idx == 4:
                    reason = "Is achiral because it is a planar molecule with a center of inversion."

        if is_active:
            code_active_indices.append(idx)
        reasons[idx] = reason

    # Compare the code's findings with the LLM's answer
    if sorted(llm_active_indices) == sorted(code_active_indices):
        return "Correct"
    else:
        error_message = "The answer is incorrect.\n"
        error_message += f"LLM's answer implies compounds {sorted(llm_active_indices)} are optically active.\n"
        error_message += f"Code analysis indicates compounds {sorted(code_active_indices)} are optically active.\n\n"
        
        # Find discrepancies
        llm_set = set(llm_active_indices)
        code_set = set(code_active_indices)
        
        wrongly_included = llm_set - code_set
        missed_by_llm = code_set - llm_set

        if wrongly_included:
            for idx in wrongly_included:
                error_message += f"Constraint violated: Compound {idx} ({compounds[idx]['name']}) was incorrectly identified as optically active. Reason: {reasons[idx]}\n"
        
        if missed_by_llm:
            for idx in missed_by_llm:
                error_message += f"Constraint violated: Compound {idx} ({compounds[idx]['name']}) was missed. Reason: {reasons[idx]}\n"
        
        return error_message.strip()

# Execute the check and print the result
result = check_answer()
print(result)