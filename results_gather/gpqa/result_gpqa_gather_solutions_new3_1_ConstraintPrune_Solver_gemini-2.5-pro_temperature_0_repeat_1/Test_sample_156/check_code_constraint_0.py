import re

def check_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function evaluates the answer based on key scientific principles:
    1.  The nature of a retrovirus (RNA genome).
    2.  The definition of a "molecular" diagnostic kit (detects nucleic acids).
    3.  The requirements for a "quick and accurate" test in an outbreak scenario.
    """
    
    # Extract the final letter answer from the text
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect", "The answer format is wrong. It should end with '<<<X>>>' where X is A, B, C, or D."
    
    answer = match.group(1)

    # Define the properties of each option based on the question
    options = {
        'A': {
            'identification': 'DNA sequencing',
            'kit': 'PCR',
            'is_molecular': True,
            'is_id_correct_for_retrovirus': False, # Incorrect: retrovirus is RNA, not DNA
            'is_kit_quick_for_outbreak': True,
            'reason': "The identification method 'DNA sequencing' is incorrect for a retrovirus, which has an RNA genome. The viral RNA must first be converted to cDNA."
        },
        'B': {
            'identification': 'cDNA sequencing',
            'kit': 'real time PCR',
            'is_molecular': True,
            'is_id_correct_for_retrovirus': True, # Correct: RNA -> cDNA -> sequence
            'is_kit_quick_for_outbreak': True, # Real-time PCR is the gold standard for speed and accuracy
            'reason': None
        },
        'C': {
            'identification': 'IgG antibodies',
            'kit': 'ELISA',
            'is_molecular': False, # Incorrect: ELISA is a serological/immunological test
            'is_id_correct_for_retrovirus': False, # Does not identify the virus's genetic material
            'is_kit_quick_for_outbreak': False, # Incorrect: IgG appears late, not suitable for quick/early diagnosis
            'reason': "This describes a serological test (ELISA), not a molecular one. Furthermore, IgG antibodies appear late in an infection, making the test unsuitable for 'quick' diagnosis of an early, active infection."
        },
        'D': {
            'identification': 'symptoms',
            'kit': 'nested PCR',
            'is_molecular': True,
            'is_id_correct_for_retrovirus': False, # Incorrect: Symptoms are non-specific and cannot be used to determine a genetic sequence
            'is_kit_quick_for_outbreak': True, # Nested PCR is sensitive but generally slower and more prone to contamination than real-time PCR
            'reason': "The identification method is scientifically invalid. A virus's genetic sequence cannot be determined from patient symptoms, which is a necessary step to design a PCR kit."
        }
    }

    selected_option = options.get(answer)

    if not selected_option:
        return "Incorrect", f"The answer '{answer}' is not a valid option (A, B, C, or D)."

    # Check against the constraints of the question
    error_messages = []
    if not selected_option['is_molecular']:
        error_messages.append(f"The kit type '{selected_option['kit']}' is not a 'molecular' diagnostic as required by the question. {selected_option['reason']}")
    
    if not selected_option['is_id_correct_for_retrovirus']:
        error_messages.append(f"The identification method '{selected_option['identification']}' is not scientifically sound for a retrovirus. {selected_option['reason']}")

    if not selected_option['is_kit_quick_for_outbreak'] and 'quick' in selected_option['reason']:
         error_messages.append(f"The proposed kit is not suitable for 'quick' diagnosis during an outbreak. {selected_option['reason']}")

    if not error_messages:
        return "Correct", ""
    else:
        # If the primary reason in the dictionary covers the errors, use it. Otherwise, join the messages.
        final_reason = selected_option['reason'] if selected_option['reason'] else " ".join(error_messages)
        return "Incorrect", f"The answer {answer} is incorrect. Reason: {final_reason}"


# The final answer provided by the LLM
llm_final_answer = """
Here is a step-by-step analysis to determine the correct answer.

**1. Analyze the Core Requirements of the Task**

*   **Pathogen:** The target is a **retrovirus**. This is the most critical piece of information. Retroviruses are RNA viruses that use an enzyme called reverse transcriptase to create a DNA copy (cDNA) of their genome.
*   **Kit Type:** The request is for a **molecular diagnostic kit**. This means the test must detect the virus's genetic material (nucleic acids: RNA or DNA), not the host's immune response (like antibodies) or viral proteins.
*   **Performance:** The kit must be for **quick and accurate** diagnosis, which is essential for managing an active outbreak.

**2. Evaluate Each Option Against the Requirements**

Let's break down each choice based on these scientific principles.

*   **A) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Careful Point:** This is incorrect from the first step. A retrovirus particle contains an RNA genome. Attempting to perform direct DNA sequencing on the virus would fail. This option misses the crucial step of converting RNA to cDNA.

*   **B) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Careful Point 1 (Identification):** This is the correct initial step. To sequence an RNA virus, its RNA must first be converted into complementary DNA (cDNA) via reverse transcription. Sequencing this cDNA allows for the precise identification of the virus's genetic code.
    *   **Careful Point 2 (Detection Kit):** A real-time PCR kit (specifically, a real-time RT-PCR kit that includes the reverse transcription step) is the modern gold standard for this purpose. It is highly sensitive, specific, and provides results much faster than traditional methods, fulfilling the "quick and accurate" requirement.

*   **C) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Careful Point 1 (Kit Type):** This describes a serological/immunological test, not a *molecular* one. It detects the host's immune response, not the virus's genetic material.
    *   **Careful Point 2 (Timing):** IgG antibodies are part of the long-term immune response and can take weeks to become detectable. A test for IgG is therefore unsuitable for *quick* diagnosis of an early, active infection. It is better for determining past exposure.

*   **D) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Careful Point:** This approach is scientifically invalid. Symptoms like fever or cough are non-specific and can be caused by many different pathogens. It is impossible to determine the unique genetic sequence of a virus—which is required to design PCR primers—based on symptoms alone.

**Conclusion**

Option B is the only one that describes a scientifically sound, logical, and modern workflow for developing a molecular diagnostic kit for a new retrovirus outbreak. It correctly accounts for the RNA nature of the virus and proposes the most appropriate technology for rapid and accurate detection.

<<<B>>>
"""

# Run the check
result, reason = check_answer(llm_final_answer)

if result == "Correct":
    print("Correct")
else:
    print(reason)
