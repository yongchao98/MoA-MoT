def check_immunology_question():
    """
    This function checks the correctness of the answer to the immunology question
    by codifying the key concepts and comparing them.
    """
    # Key information extracted from the question
    question_context = {
        "location": "secondary lymphoid organ",  # Peyer's patches
        "trigger": "after antigen encounter",
        "observation": "high variability in variable heavy chain gene"
    }

    # Database of immunological processes and their characteristics
    processes = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",  # e.g., bone marrow
            "trigger": "before antigen encounter",
            "genetic_effect": "creates initial diversity in variable region"
        },
        "B": {
            "name": "somatic hypermutation",
            "location": "secondary lymphoid organ",  # e.g., germinal centers
            "trigger": "after antigen encounter",
            "genetic_effect": "introduces high variability in variable heavy chain gene"
        },
        "C": {
            "name": "complement activation",
            "location": "serum / extracellular",
            "trigger": "after antigen encounter",
            "genetic_effect": "none (protein cascade, not genetic modification of lymphocytes)"
        },
        "D": {
            "name": "class switching recombination",
            "location": "secondary lymphoid organ",
            "trigger": "after antigen encounter",
            "genetic_effect": "changes the constant heavy chain gene, not the variable region"
        }
    }

    # The proposed answer from the LLM
    llm_answer = "B"

    # Retrieve the characteristics of the process corresponding to the answer
    answer_char = processes[llm_answer]

    # Check if the characteristics match the question's context
    if answer_char["location"] != question_context["location"]:
        return (f"Incorrect. The process in the question occurs in a {question_context['location']} (Peyer's patch), "
                f"whereas {answer_char['name']} occurs in a {answer_char['location']}.")

    if answer_char["trigger"] != question_context["trigger"]:
        return (f"Incorrect. The process in the question occurs {question_context['trigger']}, "
                f"whereas {answer_char['name']} occurs {answer_char['trigger']}.")

    # This is the most critical check
    if answer_char["genetic_effect"] != question_context["observation"]:
        return (f"Incorrect. The question observes '{question_context['observation']}'. "
                f"However, {answer_char['name']} {answer_char['genetic_effect']}.")

    # If all conditions are met, the answer is correct.
    return "Correct"

# Run the checker
result = check_immunology_question()
print(result)