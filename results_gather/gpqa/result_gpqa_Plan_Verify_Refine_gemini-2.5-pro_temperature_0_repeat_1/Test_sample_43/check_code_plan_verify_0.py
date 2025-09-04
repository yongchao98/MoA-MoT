def check_mitochondrial_experiment_answer(llm_answer: str):
    """
    Checks the correctness of the answer for the mitochondrial experiment question.

    The function evaluates four experimental options to determine which one is NOT suitable
    for studying a drug's effect on mitochondria.

    Args:
        llm_answer: The answer provided by the LLM, expected in the format "<<<X>>>".

    Returns:
        A string indicating "Correct" or an explanation of why the answer is incorrect.
    """

    # Define the biological facts for each experimental option.
    # The key is the option letter. The value is a dictionary containing its relevance
    # to mitochondrial studies and the scientific reason.
    experimental_facts = {
        'A': {
            "is_relevant": True,
            "reason": "Mito-RTP is a probe for mitochondrial reactive oxygen species (ROS). Measuring ROS is a valid way to assess mitochondrial function and stress."
        },
        'B': {
            "is_relevant": False,
            "reason": "Mitochondria do not directly take up glucose. Glucose is converted to pyruvate in the cytoplasm during glycolysis, and it is pyruvate that enters the mitochondria. Therefore, a glucose uptake assay on isolated mitochondria is a flawed experimental design."
        },
        'C': {
            "is_relevant": True,
            "reason": "The luciferase/luciferin reaction measures ATP levels. Since mitochondria are the primary producers of cellular ATP, this is a valid assay to assess mitochondrial energy production."
        },
        'D': {
            "is_relevant": True,
            "reason": "The chemical 5,5',6,6'-Tetrachloro-1,1',3,3'-tetraethylbenzimidazolylcarbocyanine iodide is JC-1, a dye used to measure mitochondrial membrane potential, which is a key indicator of mitochondrial health."
        }
    }

    # The question asks for the experiment that will NOT help.
    # We identify the correct option based on our facts.
    correct_option = None
    for option, details in experimental_facts.items():
        if not details["is_relevant"]:
            correct_option = option
            break

    # Extract the letter from the LLM's answer string (e.g., "<<<B>>>" -> "B").
    try:
        user_answer = llm_answer.strip().replace("<<<", "").replace(">>>", "")
        if user_answer not in experimental_facts:
            raise ValueError("Invalid option")
    except (AttributeError, ValueError):
        return f"Invalid answer format. The answer must be one of {list(experimental_facts.keys())} in the format <<<X>>>."

    # Compare the user's answer with the correct option.
    if user_answer == correct_option:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed explanation.
        reason_for_user_answer_incorrectness = experimental_facts[user_answer]["reason"]
        reason_for_correct_answer = experimental_facts[correct_option]["reason"]
        
        return (
            f"The provided answer '{user_answer}' is incorrect. "
            f"The experiment in option {user_answer} is a valid method for studying mitochondrial function. "
            f"Reason: {reason_for_user_answer_incorrectness}\n\n"
            f"The correct answer is '{correct_option}' because this experiment is not suitable for studying mitochondria. "
            f"Reason: {reason_for_correct_answer}"
        )

# The answer provided by the other LLM
llm_response = "<<<B>>>"

# Check the correctness of the answer
result = check_mitochondrial_experiment_answer(llm_response)
print(result)