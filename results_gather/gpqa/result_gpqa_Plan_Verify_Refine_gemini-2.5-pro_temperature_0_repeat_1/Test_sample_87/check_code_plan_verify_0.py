import collections

def check_answer():
    """
    Checks the correctness of the answer to the immunology question.

    The function models the key characteristics of the immunological processes
    listed in the options and compares them against the constraints given in the
    question's experimental setup.
    """

    # The answer provided by the LLM
    llm_answer = "D"

    # Define the knowledge base for the immunological processes
    # Each process is defined by its location, timing, and genetic target.
    knowledge_base = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",  # Bone marrow for B cells
            "timing": "pre-antigen encounter",    # During initial lymphocyte development
            "target_gene_region": "variable"      # Assembles the variable region
        },
        "B": {
            "name": "class switching recombination",
            "location": "secondary lymphoid organ", # Germinal centers (e.g., Peyer's patches)
            "timing": "post-antigen encounter",   # In activated, proliferating B cells
            "target_gene_region": "constant"      # Changes the constant region (isotype)
        },
        "C": {
            "name": "complement activation",
            "location": "blood/tissues",          # A humoral (fluid-phase) system
            "timing": "post-antigen encounter",   # Part of innate and adaptive response
            "target_gene_region": "not applicable" # Does not involve gene recombination
        },
        "D": {
            "name": "somatic hypermutation",
            "location": "secondary lymphoid organ", # Germinal centers (e.g., Peyer's patches)
            "timing": "post-antigen encounter",   # In activated, proliferating B cells
            "target_gene_region": "variable"      # Introduces point mutations into the variable region
        }
    }

    # Extract constraints from the question's text
    question_constraints = {
        "location": "secondary lymphoid organ",  # "Peyer patches" are secondary lymphoid organs.
        "timing": "post-antigen encounter",   # "proliferating cell population" after antigen delivery implies activation.
        "target_gene_region": "variable"      # "sequence their variable heavy chain gene" shows "high variability".
    }

    # Find the option that correctly matches all constraints
    correct_option = None
    mismatches = collections.defaultdict(list)

    for option, properties in knowledge_base.items():
        is_match = True
        # Check location constraint
        if properties["location"] != question_constraints["location"]:
            is_match = False
            mismatches[option].append(f"it occurs in a {properties['location']}, not a {question_constraints['location']} like a Peyer's patch")
        
        # Check timing constraint
        if properties["timing"] != question_constraints["timing"]:
            is_match = False
            mismatches[option].append(f"it occurs {properties['timing']}, not post-antigen encounter in proliferating cells")

        # Check gene region constraint
        if properties["target_gene_region"] != question_constraints["target_gene_region"]:
            is_match = False
            mismatches[option].append(f"it affects the '{properties['target_gene_region']}' gene region, not the '{question_constraints['target_gene_region']}' region")

        if is_match:
            correct_option = option

    # Final validation
    if llm_answer == correct_option:
        return "Correct"
    else:
        if llm_answer in mismatches:
            reason = f"The answer '{llm_answer}' ({knowledge_base[llm_answer]['name']}) is incorrect because " + " and ".join(mismatches[llm_answer]) + "."
            reason += f" The correct process is Somatic Hypermutation (D), which matches all conditions: it occurs in secondary lymphoid organs, after antigen encounter, and introduces high variability into the variable gene region."
            return reason
        else:
            # This case handles if the LLM answer was somehow valid but not the best fit, or if it's not in the options.
            return f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_option}'."

# Execute the check and print the result
result = check_answer()
print(result)