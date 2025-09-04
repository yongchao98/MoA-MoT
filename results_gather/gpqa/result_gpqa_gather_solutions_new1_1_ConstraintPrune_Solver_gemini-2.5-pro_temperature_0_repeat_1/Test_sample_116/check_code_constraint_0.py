def check_correctness():
    """
    This function checks the correctness of the LLM's answer by analyzing the chirality of each molecule.
    It uses the RDKit library for programmatic chemical analysis.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.rdmolops import FindPotentialStereo
    except ImportError:
        return "RDKit library not found. Please install it (`pip install rdkit-pypi`) to run this check."

    # The final answer provided by the LLM to be checked
    llm_answer = "C"

    # Define the compounds using their SMILES representation
    compounds_smiles = {
        1: "COC(=O)c1cccc(c1[N+](=O)[O-])c2c(cccc2[N+](=O)[O-])C(=O)OC",  # dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate
        2: "CC(O)C(=O)OC",  # methyl 2-hydroxypropanoate
        3: "O=C(c1ccccc1)c2ccccc2",  # benzophenone
        4: "COC(=O)/C=C/C(=O)OC"   # dimethyl fumarate
    }
    
    compound_names = {
        1: "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
        2: "methyl 2-hydroxypropanoate",
        3: "benzophenone",
        4: "dimethyl fumarate"
    }

    def is_optically_active(smiles: str) -> tuple[bool, str]:
        """
        Analyzes a molecule's SMILES string to determine if it's optically active.
        Returns a tuple (is_active, reason).
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False, "Invalid SMILES string."

        # 1. Check for chiral centers (e.g., asymmetric carbons)
        # The includeUnassigned=True flag is important for molecules without defined stereochemistry
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if chiral_centers:
            # Note: This doesn't check for meso compounds, but for this problem, it's sufficient.
            # Methyl 2-hydroxypropanoate has one chiral center and is not meso.
            return True, f"Contains {len(chiral_centers)} chiral center(s)."

        # 2. Check for axial chirality (atropisomerism)
        potential_stereo = FindPotentialStereo(mol)
        for stereo_info in potential_stereo:
            if stereo_info.type == Chem.StereoType.Bond_Atropisomer:
                return True, "Exhibits atropisomerism (axial chirality)."

        # 3. If no sources of chirality are found, the molecule is achiral.
        # This is determined by checking for symmetry elements, which RDKit implicitly does
        # by not flagging any chirality.
        return False, "Is achiral (possesses elements of symmetry like a plane or center of inversion)."

    # Determine the correct set of optically active compounds based on analysis
    correct_active_compounds = set()
    analysis_details = {}
    for num, smiles in compounds_smiles.items():
        is_active, reason = is_optically_active(smiles)
        analysis_details[num] = f"Compound {num} ({compound_names[num]}): {reason} -> Optically Active: {is_active}"
        if is_active:
            correct_active_compounds.add(num)

    # Define the options from the question
    options = {
        "A": {1, 2, 4},
        "B": {3, 4},
        "C": {1, 2},
        "D": {2, 3}
    }

    # Get the set of compounds corresponding to the LLM's answer
    llm_answer_set = options.get(llm_answer)

    if llm_answer_set is None:
        return f"Invalid Answer Format: The provided answer '{llm_answer}' is not a valid option."

    # Compare the LLM's answer with the correct answer
    if llm_answer_set == correct_active_compounds:
        return "Correct"
    else:
        error_message = f"Incorrect. The provided answer is '{llm_answer}', which corresponds to compounds {sorted(list(llm_answer_set))}.\n"
        error_message += f"The correct set of optically active compounds is {sorted(list(correct_active_compounds))}.\n\n"
        error_message += "Reasoning based on molecular analysis:\n"
        for num in sorted(analysis_details.keys()):
            error_message += f"- {analysis_details[num]}\n"
        
        missing = correct_active_compounds - llm_answer_set
        if missing:
            error_message += f"\nThe answer failed to identify optically active compound(s): {sorted(list(missing))}.\n"
        
        extra = llm_answer_set - correct_active_compounds
        if extra:
            error_message += f"The answer incorrectly identified optically inactive compound(s) as active: {sorted(list(extra))}.\n"
            
        return error_message

# Run the check and print the result
print(check_correctness())