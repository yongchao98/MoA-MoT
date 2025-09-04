def check_immunology_answer():
    """
    Checks the correctness of the answer to the immunology question.
    """
    # Define the facts/constraints presented in the question
    question_facts = {
        "location": "secondary lymphoid organ",  # Peyer's patches
        "timing": "after antigen encounter",  # Proliferating cells post-gavage
        "gene_region_affected": "variable",  # "variable heavy chain gene"
        "effect": "high variability"  # "high variability" observed
    }

    # Define the properties of each possible process
    processes = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",
            "timing": "before antigen encounter",
            "gene_region_affected": "variable",
            "effect": "initial diversity"
        },
        "B": {
            "name": "class switching recombination",
            "location": "secondary lymphoid organ",
            "timing": "after antigen encounter",
            "gene_region_affected": "constant",
            "effect": "change antibody isotype"
        },
        "C": {
            "name": "complement activation",
            "location": "blood/tissues",
            "timing": "during infection",
            "gene_region_affected": "none",
            "effect": "protein cascade for pathogen clearance"
        },
        "D": {
            "name": "somatic hypermutation",
            "location": "secondary lymphoid organ",
            "timing": "after antigen encounter",
            "gene_region_affected": "variable",
            "effect": "high variability" # for affinity maturation
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer = "D"

    # Get the properties of the process corresponding to the LLM's answer
    chosen_process = processes.get(llm_answer)

    if not chosen_process:
        return f"Invalid answer option: '{llm_answer}'. The options are A, B, C, D."

    # Check each fact against the chosen process's properties
    if chosen_process["location"] != question_facts["location"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because it occurs in a "
                f"'{chosen_process['location']}', but the question describes an event in a "
                f"'{question_facts['location']}' (Peyer's patches).")

    if chosen_process["timing"] != question_facts["timing"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because it occurs "
                f"'{chosen_process['timing']}', but the question describes a process that happens "
                f"'{question_facts['timing']}'.")

    if chosen_process["gene_region_affected"] != question_facts["gene_region_affected"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong because it affects the "
                f"'{chosen_process['gene_region_affected']}' region, but the question specifies "
                f"high variability in the '{question_facts['gene_region_affected']}' region.")

    # The 'effect' is a direct match in this case.
    # This check confirms the core observation.
    if chosen_process["effect"] != question_facts["effect"]:
        # This is a softer check, but still useful
        pass # In this specific case, the terms are identical.

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_immunology_answer()
print(result)