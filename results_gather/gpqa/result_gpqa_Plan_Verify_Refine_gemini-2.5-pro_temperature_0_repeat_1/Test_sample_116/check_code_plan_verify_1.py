def check_answer_correctness():
    """
    This function checks the correctness of the LLM's answer regarding optical isomerism.
    It uses the RDKit cheminformatics library to analyze the chirality of each molecule.

    Note: This code requires the RDKit library. You can install it via pip:
    pip install rdkit-pypi
    """
    try:
        from rdkit import Chem
    except ImportError:
        return "Error: RDKit library is required for this check but is not installed. Please run 'pip install rdkit-pypi'."

    # The LLM's answer states that compounds 1 and 2 are optically active, and 3 and 4 are not.
    # This corresponds to option B.
    llm_conclusion = {
        1: True,   # Optically active
        2: True,   # Optically active
        3: False,  # Not optically active
        4: False   # Not optically active
    }
    llm_final_option = "B"

    # Define molecules using SMILES strings, a standard chemical representation.
    smiles_dict = {
        1: "O=C(OC)c1c(cccc1[N+](=O)[O-])c2c(C(=O)OC)cccc2[N+](=O)[O-]",  # dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate
        2: "CC(O)C(=O)OC",                                              # methyl 2-hydroxypropanoate
        3: "O=C(c1ccccc1)c2ccccc2",                                     # benzophenone
        4: "COC(=O)/C=C/C(=O)OC"                                        # dimethyl fumarate
    }

    # Dictionary to store the results of our analysis.
    analysis_results = {}
    error_log = []

    for idx, smiles in smiles_dict.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            error_log.append(f"Could not parse SMILES for compound {idx}.")
            continue

        is_chiral = False

        # Check 1: Presence of tetrahedral chiral centers (e.g., a carbon with 4 different substituents).
        # This is the most common source of chirality.
        # `includeUnassigned=True` finds potential centers even if stereochemistry isn't specified.
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(chiral_centers) > 0:
            is_chiral = True

        # Check 2: Presence of axial chirality (atropisomerism).
        # This occurs in molecules like substituted biphenyls where rotation around a single bond is hindered.
        # RDKit's `FindPotentialStereo` can identify bonds that could lead to atropisomerism.
        if not is_chiral:
            stereo_info = Chem.FindPotentialStereo(mol)
            for si in stereo_info:
                # We are specifically looking for atropisomerism, which is relevant for compound 1.
                if si.type == Chem.StereoType.Bond_Atropisomer:
                    is_chiral = True
                    break
        
        analysis_results[idx] = is_chiral

    # Compare our analysis with the LLM's conclusion.
    if error_log:
        return f"An error occurred during analysis: {'; '.join(error_log)}"

    mismatches = []
    for i in range(1, 5):
        if analysis_results.get(i) != llm_conclusion.get(i):
            mismatches.append(i)

    if not mismatches:
        # The reasoning is correct. Now check if the final option matches.
        # Option B is "1 and 2".
        active_compounds = {k for k, v in analysis_results.items() if v}
        if active_compounds == {1, 2} and llm_final_option == "B":
            return "Correct"
        else:
            return f"The reasoning about which compounds are active is correct ({active_compounds}), but the final selected option '{llm_final_option}' is inconsistent."
    else:
        reason = "The answer is incorrect because its analysis of optical activity is flawed.\n"
        for idx in mismatches:
            llm_claim = 'optically active' if llm_conclusion[idx] else 'not optically active'
            our_finding = 'optically active' if analysis_results[idx] else 'not optically active'
            reason += f"For compound {idx}, the answer claims it is {llm_claim}, but chemical analysis shows it is {our_finding}.\n"
        return reason

# Execute the check and print the result.
result = check_answer_correctness()
print(result)