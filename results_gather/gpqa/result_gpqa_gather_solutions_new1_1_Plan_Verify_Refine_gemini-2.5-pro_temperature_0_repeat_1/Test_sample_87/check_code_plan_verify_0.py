def check_answer_correctness():
    """
    Checks the correctness of the answer to the immunology question.
    
    The function models the key facts from the question as constraints and
    compares them against a database of the biological processes listed
    in the options.
    """
    
    # The final answer provided by the LLM analysis
    final_answer = "B"

    # 1. Define the constraints derived from the question's text
    question_constraints = {
        "location": "secondary_lymphoid_organ",  # "Peyer patches"
        "timing": "post_antigen_encounter",      # Response to delivered proteins
        "genetic_target": "variable_region",     # "high variability" in "variable heavy chain gene"
        "cell_type": "B_cell"                    # "variable heavy chain gene" is from B cells
    }

    # 2. Define the properties of each biological process in the options
    processes = {
        "A": {
            "name": "complement activation",
            "location": "blood/tissue",
            "timing": "post_antigen_encounter",
            "genetic_target": "none",
            "cell_type": "none" # It's a protein cascade
        },
        "B": {
            "name": "somatic hypermutation",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "genetic_target": "variable_region",
            "cell_type": "B_cell"
        },
        "C": {
            "name": "VDJ recombination",
            "location": "primary_lymphoid_organ", # Occurs in bone marrow
            "timing": "pre_antigen_encounter",
            "genetic_target": "variable_region",
            "cell_type": "B_cell"
        },
        "D": {
            "name": "class switching recombination",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "genetic_target": "constant_region", # Affects the constant, not variable, region
            "cell_type": "B_cell"
        }
    }

    # 3. Check if the provided answer satisfies all constraints
    if final_answer not in processes:
        return f"Incorrect. The answer '{final_answer}' is not one of the valid options (A, B, C, D)."

    chosen_process = processes[final_answer]

    # Check location constraint
    if chosen_process["location"] != question_constraints["location"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong. "
                f"The question describes a process in a secondary lymphoid organ (Peyer's patches), "
                f"but {chosen_process['name']} occurs in a {chosen_process['location']}.")

    # Check timing constraint
    if chosen_process["timing"] != question_constraints["timing"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong. "
                f"The question describes a response after antigen encounter, "
                f"but {chosen_process['name']} occurs {chosen_process['timing']}.")

    # Check genetic target constraint
    if chosen_process["genetic_target"] != question_constraints["genetic_target"]:
        return (f"Incorrect. The answer '{chosen_process['name']}' is wrong. "
                f"The question specifies high variability in the 'variable region', "
                f"but {chosen_process['name']} primarily affects the '{chosen_process['genetic_target']}'.")

    # If all checks pass, the answer is correct
    return "Correct"

# Execute the check and print the result
result = check_answer_correctness()
print(result)