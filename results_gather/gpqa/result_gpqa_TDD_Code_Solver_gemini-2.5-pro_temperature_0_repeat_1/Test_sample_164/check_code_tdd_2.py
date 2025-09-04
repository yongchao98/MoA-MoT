import re

def check_correctness_of_llm_answer(llm_answer_string: str) -> str:
    """
    Checks the correctness of the LLM's answer for the ethylene polymerization question.

    The function establishes the ground truth for each statement based on chemical principles
    and industrial practices, then compares it with the LLM's provided answer.
    """

    # Ground truth based on chemical literature and industrial practice.
    # The question describes a tandem catalysis system for producing LLDPE from ethylene.
    ground_truth = {
        'A': {
            'is_correct': True,
            'reason': "Statement A is correct. Group VIa (Group 6) metals, especially chromium, are used in catalysts that selectively oligomerize ethylene to alpha-olefins (e.g., 1-hexene), which form the branches."
        },
        'B': {
            'is_correct': True,
            'reason': "Statement B is correct. This technology is commercialized in processes like Dow's UNIPOL II, which uses a tandem catalyst system in a single reactor."
        },
        'C': {
            'is_correct': True,
            'reason': "Statement C is correct. Catalysts based on noble metals (e.g., Palladium) for olefin polymerization exist, but their high cost is a major barrier for bulk polymer production, making this a factually sound statement."
        },
        'D': {
            'is_correct': False,
            'reason': "Statement D is incorrect. Aluminum-based activators (e.g., MAO, alkylaluminums) are widely used and often essential for activating both polymerization and oligomerization catalysts. The claim that they 'do not work' is false."
        }
    }

    # Extract the list of chosen options from the LLM's answer string.
    try:
        match = re.search(r'<<<(.+?)>>>', llm_answer_string, re.DOTALL)
        if not match:
            return "Parsing Error: The answer is not in the expected format '<<<...>>>'."
        
        # Use eval in a controlled manner to parse the list string
        parsed_answer = eval(match.group(1))
        
        if not isinstance(parsed_answer, list):
            return f"Parsing Error: The extracted answer is not a list, but a {type(parsed_answer)}."
            
        llm_selected_options = set(parsed_answer)
    except Exception as e:
        return f"Parsing Error: Failed to evaluate the answer string. Error: {e}"

    # Determine the set of correct options from the ground truth.
    true_options = {key for key, value in ground_truth.items() if value['is_correct']}

    # Compare the LLM's answer with the ground truth.
    if llm_selected_options == true_options:
        return "Correct"
    else:
        # Generate a detailed report of the discrepancy.
        error_messages = []
        
        # Check for correct options that the LLM missed.
        missing_options = true_options - llm_selected_options
        for option in sorted(list(missing_options)):
            error_messages.append(f"The answer is missing the correct statement '{option}'. {ground_truth[option]['reason']}")

        # Check for incorrect options that the LLM included.
        incorrectly_included_options = llm_selected_options - true_options
        for option in sorted(list(incorrectly_included_options)):
            error_messages.append(f"The answer incorrectly includes statement '{option}'. {ground_truth[option]['reason']}")
        
        # Check for any invalid options (not A, B, C, or D).
        all_valid_options = set(ground_truth.keys())
        invalid_options = llm_selected_options - all_valid_options
        for option in invalid_options:
            error_messages.append(f"The answer includes an invalid option '{option}'.")

        return "\n".join(error_messages)

# The provided response from the other LLM
llm_response = """
All tests passed
Answer: The question is ambiguous as multiple statements are correct. The set of factually correct statements is ['A', 'B', 'C'].

<<<['A', 'B', 'C']>>>
"""

# Run the check
result = check_correctness_of_llm_answer(llm_response)
print(result)