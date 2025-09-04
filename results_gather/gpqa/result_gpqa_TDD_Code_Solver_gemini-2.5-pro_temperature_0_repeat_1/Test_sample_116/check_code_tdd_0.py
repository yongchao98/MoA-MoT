import re

def check_answer():
    """
    Checks the correctness of the given answer for the organic chemistry question.
    """
    # Correct analysis of each compound's optical activity based on chemical principles.
    # True = Optically Active, False = Optically Inactive
    correct_analysis = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_active": True,
            "reason": "Exhibits atropisomerism due to restricted rotation around the biphenyl C-C bond caused by bulky ortho substituents, making it chiral."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_active": True,
            "reason": "Contains a chiral center (the carbon bonded to -H, -OH, -CH3, and -COOCH3)."
        },
        3: {
            "name": "benzophenone",
            "is_active": False,
            "reason": "Is achiral because it possesses a plane of symmetry."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_active": False,
            "reason": "Is achiral because it is planar and has a center of inversion."
        }
    }

    # Determine the set of optically active compounds from the correct analysis
    active_compounds_indices = {key for key, value in correct_analysis.items() if value["is_active"]}

    # Define the options given in the question
    options = {
        "A": {1, 2},
        "B": {2, 3},
        "C": {1, 2, 4},
        "D": {3, 4}
    }

    # Find the correct option letter
    correct_option = None
    for option, indices in options.items():
        if indices == active_compounds_indices:
            correct_option = option
            break

    # The provided answer from the LLM
    llm_answer_text = "<<<A>>>"
    
    # Extract the letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D, but got '{llm_answer_text}'."
        
    llm_option = match.group(1)

    # Check if the LLM's answer matches the correct option
    if llm_option == correct_option:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrectness
        llm_selected_compounds = options.get(llm_option, set())
        
        error_messages = []
        # Check for compounds incorrectly included by the LLM
        for index in llm_selected_compounds:
            if not correct_analysis[index]["is_active"]:
                error_messages.append(f"Compound {index} ({correct_analysis[index]['name']}) was incorrectly identified as optically active. Reason: {correct_analysis[index]['reason']}")

        # Check for compounds incorrectly excluded by the LLM
        for index in active_compounds_indices:
            if index not in llm_selected_compounds:
                error_messages.append(f"Compound {index} ({correct_analysis[index]['name']}) was incorrectly identified as not optically active. Reason: {correct_analysis[index]['reason']}")

        reason = f"The answer '{llm_option}' is incorrect. The correct answer is '{correct_option}'.\n"
        reason += "The set of optically active compounds should be {1, 2}.\n"
        reason += "\n".join(error_messages)
        return reason

# Run the checker
result = check_answer()
print(result)