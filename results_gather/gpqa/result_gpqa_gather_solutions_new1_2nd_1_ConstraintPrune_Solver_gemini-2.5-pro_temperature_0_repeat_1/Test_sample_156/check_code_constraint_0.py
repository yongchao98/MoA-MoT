def check_diagnostic_kit_design():
    """
    Checks the correctness of the selected answer for designing a retrovirus diagnostic kit.

    The function evaluates each option based on three key constraints from the question:
    1.  Kit Type: Must be "molecular" (detects nucleic acids).
    2.  Pathogen: Must be for a "retrovirus" (RNA genome, requires reverse transcription).
    3.  Goal: Must be "quick and accurate" (favors modern, rapid, sensitive methods).
    """

    # The options as presented in the question prompt.
    options = {
        'A': "First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis",
        'B': "First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.",
        'C': "First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.",
        'D': "First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis."
    }

    # The final answer to be checked.
    final_answer = 'B'

    # Analysis of each option against the constraints.
    analysis = {
        'A': {
            'is_correct': False,
            'reason': "This is incorrect because a retrovirus has an RNA genome. Direct 'DNA sequencing' would fail. It misunderstands the fundamental biology of the pathogen."
        },
        'B': {
            'is_correct': True,
            'reason': "This is the correct approach. It correctly identifies the need to convert the viral RNA to cDNA for sequencing ('cDNA sequencing') and uses the modern gold standard for quick and accurate molecular diagnosis ('real time PCR')."
        },
        'C': {
            'is_correct': False,
            'reason': "This is incorrect for two reasons: 1) It describes an ELISA kit, which is an immunological test, not a 'molecular' kit. 2) It targets IgG antibodies, which are markers for late-stage or past infection, not for 'quick' early diagnosis."
        },
        'D': {
            'is_correct': False,
            'reason': "This is incorrect because the initial step of identifying a virus based on symptoms is scientifically invalid for designing a specific molecular test. The genetic sequence is required and cannot be derived from symptoms."
        }
    }

    # Check if the selected final answer is correct.
    if final_answer in analysis:
        if analysis[final_answer]['is_correct']:
            return "Correct"
        else:
            # The provided answer is incorrect. Find the correct one and explain the error.
            correct_option = None
            for option, details in analysis.items():
                if details['is_correct']:
                    correct_option = option
                    break
            
            return (f"Incorrect. The final answer '{final_answer}' is wrong. "
                    f"Reason: {analysis[final_answer]['reason']} "
                    f"The correct answer is '{correct_option}' because it is the only option that satisfies all constraints: it is a molecular method, it correctly handles the RNA genome of a retrovirus, and it uses a technology suitable for quick and accurate diagnosis.")
    else:
        return f"Invalid option '{final_answer}' provided. The options are A, B, C, D."

# Execute the check and print the result.
result = check_diagnostic_kit_design()
print(result)