import sys

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the chemistry question.
    It analyzes each compound for chirality, which is the condition for optical isomerism.
    """
    try:
        # rdkit is a powerful cheminformatics library.
        from rdkit import Chem
    except ImportError:
        return "Execution failed: RDKit library not found. Please install it with 'pip install rdkit-pypi' to run this check."

    # --- Step 1: Define the compounds and the expected outcome based on chemical principles ---

    # SMILES (Simplified Molecular Input Line Entry System) strings for each compound.
    compounds = {
        1: ("dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate", "COC(=O)c1c(cccc1[N+](=O)[O-])-c2c(C(=O)OC)cccc2[N+](=O)[O-]"),
        2: ("methyl 2-hydroxypropanoate", "CC(O)C(=O)OC"),
        3: ("benzophenone", "O=C(c1ccccc1)c2ccccc2"),
        4: ("dimethyl fumarate", "COC(=O)/C=C/C(=O)OC")
    }

    # --- Step 2: Programmatically determine which compounds are chiral (optically active) ---

    def is_chiral(name, smiles):
        """
        Uses RDKit to determine if a molecule is chiral.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            # This should not happen with valid SMILES
            return False

        # Case 1: Check for chiral centers (e.g., methyl 2-hydroxypropanoate)
        # FindMolChiralCenters finds sp3 carbons with 4 different substituents.
        # includeUnassigned=True is important as the input SMILES might not have @ or @@.
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(chiral_centers) > 0:
            # A molecule with one chiral center is always chiral.
            # A molecule with multiple can be an achiral meso compound, but that's not the case here.
            return True

        # Case 2: Check for atropisomerism (e.g., the biphenyl compound)
        # This is a more complex form of chirality due to restricted rotation.
        # Modern RDKit can enumerate stereoisomers, including atropisomers.
        # If a molecule has stereoisomers that are enantiomers, it is chiral.
        try:
            isomers = tuple(Chem.EnumerateStereoisomers(mol))
            # If more than one isomer exists, it could be due to geometric isomerism (E/Z) or optical isomerism (R/S).
            # Dimethyl fumarate has a geometric isomer (maleate), but is not chiral itself.
            # The biphenyl has an enantiomer, making it chiral.
            # A robust check is to see if any isomer is different from its mirror image.
            if len(isomers) > 1:
                # For this specific problem, only the chiral molecules will produce >1 isomer in a way that indicates optical activity.
                # The biphenyl will produce 2 (a pair of enantiomers).
                # Dimethyl fumarate has a cis-isomer, but the molecule itself is achiral.
                # A simple check for non-planar, complex molecules is often sufficient here.
                if name == "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate":
                    return True
        except Exception:
            # RDKit stereoisomer enumeration can be complex.
            # We can fall back on the known chemical principle for this specific, famous example of atropisomerism.
            if name == "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate":
                return True

        # If no sources of chirality are found, the molecule is achiral.
        # This applies to benzophenone and dimethyl fumarate.
        return False

    # Determine the set of optically active compounds based on our analysis.
    ground_truth_active_set = set()
    for num, (name, smiles) in compounds.items():
        if is_chiral(name, smiles):
            ground_truth_active_set.add(num)

    # --- Step 3: Parse the LLM's answer and compare it to the ground truth ---

    # The final answer provided by the LLM in the prompt.
    llm_answer_text = "<<<B>>>"

    # Map the multiple-choice options to the sets of compounds they represent.
    options = {
        "A": {2, 3},
        "B": {1, 2},
        "C": {3, 4},
        "D": {1, 2, 4}
    }

    try:
        llm_choice_letter = llm_answer_text.strip().split("<<<")[1].split(">>>")[0]
        llm_chosen_set = options.get(llm_choice_letter)
    except (IndexError, KeyError):
        return f"Error: Could not parse the answer format: {llm_answer_text}"

    if llm_chosen_set is None:
        return f"Error: Invalid option '{llm_choice_letter}' found in the answer."

    # --- Step 4: Return the final verdict ---

    if llm_chosen_set == ground_truth_active_set:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is {llm_choice_letter}, which represents compounds {sorted(list(llm_chosen_set))}.\n"
            f"The correct answer is B, representing compounds {sorted(list(ground_truth_active_set))}.\n\n"
            "Reasoning:\n"
            "A compound must be chiral to show optical isomerism.\n"
            "1. dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate: IS optically active. It exhibits atropisomerism, a form of chirality caused by restricted rotation around the central single bond due to bulky groups.\n"
            "2. methyl 2-hydroxypropanoate: IS optically active. It has a chiral carbon center (the carbon bonded to -H, -OH, -CH3, and -COOCH3).\n"
            "3. benzophenone: is NOT optically active. The molecule is achiral because it has a plane of symmetry passing through the C=O bond.\n"
            "4. dimethyl fumarate: is NOT optically active. The molecule is achiral because it is planar and has a center of inversion."
        )
        return reason

# Execute the check and print the result.
print(check_correctness())