import re

def check_optical_isomerism_answer():
    """
    This function checks the correctness of the answer for the optical isomerism question.
    It encodes the chemical properties of each compound and compares the provided answer against the ground truth.
    """
    # Ground truth based on chemical principles
    # A compound shows optical isomerism if it is chiral.
    compounds_properties = {
        1: {
            'name': "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            'is_chiral': True,
            'reason': "It exhibits atropisomerism. The bulky ortho groups (-NO2, -COOCH3) restrict rotation around the C-C single bond, locking the molecule in a non-planar, chiral conformation."
        },
        2: {
            'name': 'methyl 2-hydroxypropanoate',
            'is_chiral': True,
            'reason': "It has a chiral center. The carbon at position 2 is bonded to four different groups: -H, -OH, -CH3, and -COOCH3."
        },
        3: {
            'name': 'benzophenone',
            'is_chiral': False,
            'reason': "It is achiral. The molecule has a plane of symmetry that passes through the C=O bond and bisects the angle between the two phenyl rings."
        },
        4: {
            'name': 'dimethyl fumarate',
            'is_chiral': False,
            'reason': "It is achiral. The molecule is planar (possessing a plane of symmetry) and has a center of inversion at the midpoint of the C=C double bond."
        }
    }

    # Define the options as sets of compound numbers
    options = {
        'A': {3, 4},
        'B': {2, 3},
        'C': {1, 2},
        'D': {1, 2, 4}
    }

    # The final answer provided by the LLM to be checked
    llm_answer_text = "<<<C>>>"

    # Extract the letter from the answer string
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format: '{llm_answer_text}'. Expected format is <<<X>>> where X is A, B, C, or D."

    llm_option = match.group(1)
    llm_answered_set = options.get(llm_option)

    # Determine the correct set of compounds that are chiral
    correct_chiral_set = {num for num, props in compounds_properties.items() if props['is_chiral']}

    # Check if the LLM's answer matches the correct set
    if llm_answered_set == correct_chiral_set:
        return "Correct"
    else:
        # If incorrect, build a detailed reason
        error_message = f"Incorrect. The provided answer is '{llm_option}', which corresponds to compounds {llm_answered_set}.\n"
        error_message += f"The correct set of optically active compounds is {correct_chiral_set}.\n\n"
        error_message += "Here is a breakdown of the errors:\n"

        # Check for compounds incorrectly included
        incorrectly_included = llm_answered_set - correct_chiral_set
        for num in sorted(list(incorrectly_included)):
            error_message += f"- Compound {num} ('{compounds_properties[num]['name']}') was incorrectly identified as optically active. It is achiral.\n"
            error_message += f"  Reason: {compounds_properties[num]['reason']}\n"

        # Check for compounds incorrectly excluded
        incorrectly_excluded = correct_chiral_set - llm_answered_set
        for num in sorted(list(incorrectly_excluded)):
            error_message += f"- Compound {num} ('{compounds_properties[num]['name']}') was incorrectly identified as not optically active. It is chiral.\n"
            error_message += f"  Reason: {compounds_properties[num]['reason']}\n"
            
        return error_message.strip()

# Execute the check and print the result
result = check_optical_isomerism_answer()
print(result)