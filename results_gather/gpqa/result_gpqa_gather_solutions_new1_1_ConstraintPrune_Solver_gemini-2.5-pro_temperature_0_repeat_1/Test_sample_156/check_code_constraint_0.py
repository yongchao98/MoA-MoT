def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer for designing a retrovirus diagnostic kit.

    The function hardcodes the question's constraints and the options' properties.
    It then evaluates which option is scientifically correct and compares it
    to the provided answer.
    """
    # The final answer provided by the LLM analysis.
    llm_answer = "A"

    # --- Step 1: Define the constraints from the question ---
    constraints = {
        "pathogen_type": "retrovirus",
        "pathogen_genome": "RNA",
        "kit_type": "molecular",  # Must detect nucleic acids (RNA/DNA)
        "goal": ["quick", "accurate"]
    }

    # --- Step 2: Define the properties of each option ---
    options_analysis = {
        "A": {
            "description": "cDNA sequencing -> real time PCR kit",
            "is_molecular": True,
            "is_correct_for_rna_virus": True, # cDNA sequencing is correct for RNA
            "is_quick_for_early_detection": True, # Real-time PCR is fast
            "is_scientifically_sound": True
        },
        "B": {
            "description": "DNA sequencing -> PCR kit",
            "is_molecular": True,
            "is_correct_for_rna_virus": False, # Direct DNA sequencing is wrong for an RNA virus
            "is_quick_for_early_detection": True,
            "is_scientifically_sound": False
        },
        "C": {
            "description": "Symptoms -> nested PCR kit",
            "is_molecular": True,
            "is_correct_for_rna_virus": False, # Cannot design a kit from symptoms
            "is_quick_for_early_detection": False, # Design is impossible/not quick
            "is_scientifically_sound": False
        },
        "D": {
            "description": "IgG antibodies -> ELISA kit",
            "is_molecular": False, # ELISA is an immunological test
            "is_correct_for_rna_virus": True, # Not applicable, but doesn't violate this
            "is_quick_for_early_detection": False, # IgG is a late-stage antibody
            "is_scientifically_sound": False # Fails on two other constraints
        }
    }

    # --- Step 3: Systematically evaluate options to find the correct one ---
    correct_option = None
    for option, analysis in options_analysis.items():
        # Check against all constraints
        if (analysis["is_molecular"] and
            analysis["is_correct_for_rna_virus"] and
            analysis["is_quick_for_early_detection"] and
            analysis["is_scientifically_sound"]):
            correct_option = option
            break

    # --- Step 4: Compare the determined correct answer with the LLM's answer ---
    if llm_answer == correct_option:
        return "Correct"
    else:
        # Provide a reason why the LLM's answer is incorrect
        chosen_analysis = options_analysis.get(llm_answer)
        if not chosen_analysis:
            return f"Incorrect: The provided answer '{llm_answer}' is not a valid option."

        reasons = []
        if not chosen_analysis["is_molecular"]:
            reasons.append("it describes an immunological test (ELISA), not a molecular one as required.")
        if not chosen_analysis["is_correct_for_rna_virus"]:
            if llm_answer == 'B':
                reasons.append("it incorrectly suggests direct DNA sequencing for a retrovirus, which has an RNA genome.")
            elif llm_answer == 'C':
                reasons.append("it suggests designing a molecular test based on symptoms, which is scientifically impossible.")
        if not chosen_analysis["is_quick_for_early_detection"]:
            reasons.append("it is not suitable for 'quick' early diagnosis because it relies on late-stage IgG antibodies.")
        
        if not reasons:
             reasons.append("it fails to meet all the scientific requirements of the task.")

        return f"Incorrect: The correct answer is {correct_option}, but the provided answer is {llm_answer}. The given answer is wrong because {' and '.join(reasons)}"

# Execute the check and print the result
result = check_answer_correctness()
print(result)