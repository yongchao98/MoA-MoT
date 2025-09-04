import collections

def check_correctness_of_answer():
    """
    This function programmatically checks the correctness of the final answer by evaluating
    each option against the scientific constraints given in the question.

    The question requires designing a "molecular diagnostic kit" for a "retrovirus"
    that allows for "quick and accurate detection".

    Constraints:
    1.  **Molecular Diagnostic**: The kit must detect nucleic acids (RNA/DNA), not proteins like antibodies.
    2.  **Retrovirus**: The pathogen has an RNA genome. Any sequencing or amplification must account for this, typically by converting RNA to complementary DNA (cDNA) first.
    3.  **Quick Detection**: The method must be suitable for early diagnosis in an outbreak, ruling out tests that detect late-stage markers.
    4.  **Valid Identification**: The initial step to identify the virus must be scientifically sound for designing a molecular test.
    """

    # The final answer provided by the LLM that we need to check.
    llm_final_answer = 'A'

    # Define the properties of each option based on the original question text.
    options_properties = {
        'A': {
            "description": "cDNA sequencing -> real time PCR kit",
            "is_molecular": True,
            "handles_retrovirus_biology": True,  # Correct: cDNA sequencing starts from RNA.
            "is_id_method_valid": True,
            "is_quick_for_outbreak": True  # Correct: Real-time PCR is fast and detects early infection.
        },
        'B': {
            "description": "DNA sequencing -> PCR kit",
            "is_molecular": True,
            "handles_retrovirus_biology": False,  # Incorrect: A retrovirus has an RNA genome.
            "is_id_method_valid": False, # The method is biologically flawed for this pathogen.
            "is_quick_for_outbreak": True
        },
        'C': {
            "description": "IgG antibodies -> ELISA kit",
            "is_molecular": False,  # Incorrect: ELISA is an immunological test, not molecular.
            "handles_retrovirus_biology": "N/A",
            "is_id_method_valid": True,
            "is_quick_for_outbreak": False  # Incorrect: IgG is a late-stage marker, not for "quick" early diagnosis.
        },
        'D': {
            "description": "Symptoms -> nested PCR kit",
            "is_molecular": True,
            "handles_retrovirus_biology": "N/A",
            "is_id_method_valid": False,  # Incorrect: Cannot design a specific molecular test from non-specific symptoms.
            "is_quick_for_outbreak": True
        }
    }

    correct_options = []
    error_messages = collections.defaultdict(list)

    for option, props in options_properties.items():
        is_valid = True
        
        # Check Constraint 1: Must be a "molecular diagnostic" kit.
        if not props["is_molecular"]:
            is_valid = False
            error_messages[option].append("it is an immunological test (ELISA), not a molecular test as required.")
            
        # Check Constraint 2: Method must be valid for a retrovirus (RNA genome).
        if not props["handles_retrovirus_biology"]:
            is_valid = False
            error_messages[option].append("it misunderstands the biology of a retrovirus, which has an RNA genome. The proposed method (direct DNA sequencing) would fail.")
            
        # Check Constraint 3: Identification method must be scientifically sound.
        if not props["is_id_method_valid"]:
            is_valid = False
            # Avoid duplicating the error from constraint 2
            if props["handles_retrovirus_biology"] is not False:
                 error_messages[option].append("the identification method (e.g., using symptoms) is scientifically invalid for designing a specific molecular test.")

        # Check Constraint 4: Must be "quick" for an outbreak (early diagnosis).
        if not props["is_quick_for_outbreak"]:
            is_valid = False
            error_messages[option].append("it is not suitable for 'quick' early diagnosis, as it detects late-stage markers (IgG).")

        if is_valid:
            correct_options.append(option)

    # Final evaluation and result generation.
    if len(correct_options) == 1:
        correct_answer = correct_options[0]
        if llm_final_answer == correct_answer:
            return "Correct"
        else:
            reason_for_llm_fail = "; ".join(error_messages[llm_final_answer])
            return (f"Incorrect. The provided answer is {llm_final_answer}, but the correct answer is {correct_answer}. "
                    f"The provided answer {llm_final_answer} is wrong because {reason_for_llm_fail}")
    elif len(correct_options) == 0:
        return "Incorrect. No option correctly satisfies all the constraints of the question."
    else:
        return f"Incorrect. The question is ambiguous as multiple options ({', '.join(correct_options)}) satisfy all constraints."

# Execute the checking function and print the result.
result = check_correctness_of_answer()
print(result)