import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the immunology question.

    The function analyzes the key facts from the question and compares them against the
    properties of the immunological processes listed in the options.
    """
    # --- Define facts from the question ---
    question_facts = {
        "location": "secondary lymphoid organ",
        "timing": "post-antigen encounter",
        "cell_state": "proliferating",
        "genetic_locus": "variable region"
    }

    # --- Define properties of the immunological processes ---
    processes = {
        "A": {
            "name": "Complement activation",
            "description": "A protein cascade, not a B cell genetic process causing variability."
        },
        "B": {
            "name": "Somatic hypermutation",
            "location": "secondary lymphoid organ",
            "timing": "post-antigen encounter",
            "cell_state": "proliferating",
            "genetic_locus": "variable region",
            "description": "Introduces point mutations in the variable region of proliferating B cells in secondary lymphoid organs after antigen encounter."
        },
        "C": {
            "name": "Class switching recombination",
            "location": "secondary lymphoid organ",
            "timing": "post-antigen encounter",
            "cell_state": "proliferating",
            "genetic_locus": "constant region",
            "description": "Affects the constant region, not the variable region as specified in the question."
        },
        "D": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",
            "timing": "pre-antigen encounter",
            "cell_state": "developing",
            "genetic_locus": "variable region",
            "description": "Occurs in primary lymphoid organs before antigen encounter, not in proliferating cells in secondary organs after."
        }
    }

    # --- Determine the correct option based on the facts ---
    correct_option = None
    for option, properties in processes.items():
        # Skip non-genetic processes
        if "location" not in properties:
            continue
        
        if all(properties.get(key) == question_facts.get(key) for key in question_facts):
            correct_option = option
            break

    # --- Validate the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer must be in the format <<<X>>> where X is A, B, C, or D."

    provided_option = match.group(1)

    if provided_option == correct_option:
        return "Correct"
    else:
        incorrect_process_info = processes.get(provided_option)
        reason = (f"The provided answer '{provided_option}' ({incorrect_process_info['name']}) is incorrect. "
                  f"{incorrect_process_info['description']}")
        return reason

# The final answer provided by the LLM to be checked
llm_final_answer = "<<<B>>>"

# Run the check
result = check_answer_correctness(llm_final_answer)
print(result)