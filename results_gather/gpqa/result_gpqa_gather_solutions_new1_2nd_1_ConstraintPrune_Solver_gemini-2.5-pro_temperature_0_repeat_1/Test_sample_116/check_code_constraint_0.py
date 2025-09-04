import re

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer by analyzing each compound for chirality.
    """

    # A simplified function to determine chirality based on established chemical rules.
    def is_chiral(compound_name):
        """
        Determines if a compound is chiral based on its name and known structural rules.
        Returns a tuple (is_chiral_boolean, reason_string).
        """
        # Rule 1: Check for atropisomerism in substituted biphenyls.
        # Biphenyls with bulky groups at all four ortho-positions (2, 2', 6, 6') are chiral.
        if "biphenyl" in compound_name and all(sub in compound_name for sub in ["6,6'-", "2,2'-"]):
            return (True, "This is a substituted biphenyl with bulky groups at all four ortho-positions, causing restricted rotation (atropisomerism). This leads to axial chirality.")

        # Rule 2: Check for a standard chiral center (a carbon with four different groups).
        if compound_name == "methyl 2-hydroxypropanoate":
            # Structure: CH3-CH(OH)-COOCH3. The carbon at position 2 is bonded to -H, -OH, -CH3, and -COOCH3.
            return (True, "This molecule contains a chiral center (the carbon at position 2), which is bonded to four different groups.")

        # Rule 3: Check for elements of symmetry that make a molecule achiral.
        if compound_name == "benzophenone":
            # Structure: (C6H5)2C=O. This molecule has a plane of symmetry.
            return (False, "This molecule is achiral because it possesses a plane of symmetry that passes through the C=O bond and bisects the two phenyl rings.")

        if compound_name == "dimethyl fumarate":
            # Structure: trans-CH3OOC-CH=CH-COOCH3. This molecule has a center of inversion.
            return (False, "This molecule is achiral. As a trans-alkene, it is planar and has a center of inversion at the midpoint of the C=C double bond.")

        # Fallback for unrecognized compounds
        return (None, "Compound structure not recognized by this checker.")

    # --- Main Logic ---
    
    # The compounds and options from the question
    compounds = {
        1: "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
        2: "methyl 2-hydroxypropanoate",
        3: "benzophenone",
        4: "dimethyl fumarate"
    }
    
    options = {
        "A": {1, 2},
        "B": {2, 3},
        "C": {1, 2, 4},
        "D": {3, 4}
    }

    # The final answer provided by the LLM to be checked
    llm_answer_text = "<<<A>>>"
    
    # Extract the letter from the answer format
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format: {llm_answer_text}"
    llm_option = match.group(1)

    # Determine the correct answer by analyzing each compound
    correctly_identified_chiral_compounds = set()
    analysis_details = []
    for num, name in compounds.items():
        chiral, reason = is_chiral(name)
        if chiral:
            correctly_identified_chiral_compounds.add(num)
        analysis_details.append(f"  - Compound {num} ({name}): {'Optically active (chiral)' if chiral else 'Not optically active (achiral)'}. Reason: {reason}")

    # Find which option corresponds to the correct set of compounds
    correct_option = None
    for option_letter, compound_set in options.items():
        if compound_set == correctly_identified_chiral_compounds:
            correct_option = option_letter
            break
            
    if correct_option is None:
        return "Checker Error: The analysis did not match any of the available options."

    # Compare the LLM's answer with the derived correct answer
    if llm_option == correct_option:
        return "Correct"
    else:
        error_message = f"Incorrect. The provided answer was {llm_option}, but the correct answer is {correct_option}.\n\n"
        error_message += "Here is the correct analysis:\n"
        error_message += "\n".join(analysis_details)
        error_message += f"\n\nConclusion: The compounds that show optical isomerism are {sorted(list(correctly_identified_chiral_compounds))}, which corresponds to option {correct_option}."
        return error_message

# Execute the check and print the result
result = check_answer_correctness()
print(result)