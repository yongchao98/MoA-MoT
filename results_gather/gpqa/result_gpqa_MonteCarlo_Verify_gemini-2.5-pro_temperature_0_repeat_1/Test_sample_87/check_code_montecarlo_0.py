def check_immunology_answer():
    """
    This function checks the correctness of the answer to the immunology question.
    It verifies the answer against a set of established facts from the question.
    """
    # --- Step 1: Define the key facts and constraints from the question ---
    # 1. Location: Peyer's patches are secondary lymphoid organs where germinal centers form.
    # 2. Cellular State: "Proliferating cell population" implies clonal expansion after antigen encounter.
    # 3. Genetic Locus: The sequencing target is the "variable heavy chain gene".
    # 4. Observation: The key finding is "high variability" in this specific gene region.
    question_constraints = {
        "location": "secondary lymphoid organ",
        "timing": "after antigen encounter",
        "genetic_locus": "variable region",
        "outcome": "high variability"
    }

    # --- Step 2: Define the properties of each biological process (the options) ---
    knowledge_base = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",  # Bone marrow for B cells
            "timing": "before antigen encounter", # During B cell development
            "genetic_locus": "variable region",
            "outcome": "initial diversity"
        },
        "B": {
            "name": "Class switching recombination",
            "location": "secondary lymphoid organ",
            "timing": "after antigen encounter",
            "genetic_locus": "constant region", # This is the key difference
            "outcome": "change of antibody isotype (e.g., IgM to IgG)"
        },
        "C": {
            "name": "Complement activation",
            "location": "N/A (protein cascade in blood/tissues)",
            "timing": "N/A (not a lymphocyte genetic process)",
            "genetic_locus": "N/A (does not alter lymphocyte genes)",
            "outcome": "pathogen lysis and opsonization"
        },
        "D": {
            "name": "Somatic hypermutation",
            "location": "secondary lymphoid organ", # Occurs in germinal centers
            "timing": "after antigen encounter",   # During B cell proliferation
            "genetic_locus": "variable region",    # Targets V genes
            "outcome": "high variability"         # Introduces point mutations, leading to affinity maturation
        }
    }

    # --- Step 3: The provided answer to be checked ---
    llm_answer = "D"

    # --- Step 4: Check the answer against the constraints ---
    selected_option = knowledge_base.get(llm_answer)

    if not selected_option:
        return f"Invalid answer option provided: {llm_answer}. Valid options are A, B, C, D."

    # Check 1: Location
    if selected_option["location"] != question_constraints["location"]:
        return (f"Incorrect. The answer '{selected_option['name']}' is wrong because it occurs in a "
                f"{selected_option['location']}, but the question describes events in Peyer's patches, "
                f"which are a {question_constraints['location']}.")

    # Check 2: Timing
    if selected_option["timing"] != question_constraints["timing"]:
        return (f"Incorrect. The answer '{selected_option['name']}' is wrong because it occurs "
                f"{selected_option['timing']}, but the question describes proliferating cells, which "
                f"indicates a process happening {question_constraints['timing']}.")

    # Check 3: Genetic Locus
    if selected_option["genetic_locus"] != question_constraints["genetic_locus"]:
        return (f"Incorrect. The answer '{selected_option['name']}' is wrong because it affects the "
                f"{selected_option['genetic_locus']}, but the question explicitly states that the "
                f"'{question_constraints['genetic_locus']}' was sequenced and showed high variability.")

    # Check 4: Outcome
    if selected_option["outcome"] != question_constraints["outcome"]:
        # This check confirms the process matches the direct observation.
        return (f"Incorrect. The answer '{selected_option['name']}' is wrong because its primary result is "
                f"'{selected_option['outcome']}', which does not directly match the key observation of "
                f"'{question_constraints['outcome']}' in the variable region.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_immunology_answer()
print(result)