import re

def check_chemistry_answer():
    """
    Checks the correctness of the answer to the chemistry question about optical isomerism.
    
    The function encodes the chemical properties of each compound to determine which ones are chiral
    and thus exhibit optical isomerism. It then compares this ground truth against the provided answer.
    """

    # --- Ground Truth Analysis ---
    # A compound shows optical isomerism if it is chiral.
    # 1. dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate: Chiral due to atropisomerism (restricted rotation).
    # 2. methyl 2-hydroxypropanoate: Chiral due to a stereocenter (C2 is bonded to 4 different groups).
    # 3. benzophenone: Achiral due to a plane of symmetry.
    # 4. dimethyl fumarate: Achiral due to a plane of symmetry and a center of inversion.
    
    compounds_properties = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_optically_active": True,
            "reason": "It exhibits atropisomerism (a form of axial chirality) due to bulky ortho-substituents that restrict rotation, making the molecule chiral."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_optically_active": True,
            "reason": "It contains a chiral center (the second carbon is bonded to four different groups: -H, -OH, -CH3, and -COOCH3)."
        },
        3: {
            "name": "benzophenone",
            "is_optically_active": False,
            "reason": "It is achiral because it possesses a plane of symmetry."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_optically_active": False,
            "reason": "It is achiral because it is planar and has a center of inversion."
        }
    }

    # --- Answer Verification ---
    
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "<<<B>>>"

    # Map options to the compound numbers they represent.
    options_map = {
        "A": {3, 4},
        "B": {1, 2},
        "C": {2, 3},
        "D": {1, 2, 4}
    }

    # Extract the letter from the answer format.
    match = re.search(r'<<<([A-D])>>>', llm_final_answer)
    if not match:
        return f"Invalid answer format: {llm_final_answer}. Expected format like '<<<A>>>'."
    
    llm_option_letter = match.group(1)
    
    if llm_option_letter not in options_map:
        return f"Invalid option '{llm_option_letter}'. Valid options are A, B, C, D."

    # Determine the set of compounds that are actually optically active based on our analysis.
    correct_active_compounds = {num for num, props in compounds_properties.items() if props["is_optically_active"]}

    # Get the set of compounds from the LLM's chosen option.
    llm_chosen_compounds = options_map[llm_option_letter]

    # Compare the correct set with the LLM's chosen set.
    if correct_active_compounds == llm_chosen_compounds:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness.
        error_messages = []
        
        # Check for compounds incorrectly included in the answer.
        incorrectly_included = llm_chosen_compounds - correct_active_compounds
        for num in sorted(list(incorrectly_included)):
            reason = compounds_properties[num]['reason']
            error_messages.append(f"Compound {num} is incorrectly included in the answer; it is not optically active because {reason}")

        # Check for compounds incorrectly excluded from the answer.
        incorrectly_excluded = correct_active_compounds - llm_chosen_compounds
        for num in sorted(list(incorrectly_excluded)):
            reason = compounds_properties[num]['reason']
            error_messages.append(f"Compound {num} is incorrectly excluded from the answer; it is optically active because {reason}")
        
        return "Incorrect. " + ". ".join(error_messages)

# Execute the check and print the result.
print(check_chemistry_answer())