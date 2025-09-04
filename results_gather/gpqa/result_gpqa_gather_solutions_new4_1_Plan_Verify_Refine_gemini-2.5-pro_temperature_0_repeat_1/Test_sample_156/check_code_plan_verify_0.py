import re

def check_answer(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function encodes the scientific principles of molecular diagnostics for retroviruses
    and evaluates the chosen option against these principles.
    """

    # --- Define Scientific Constraints ---
    constraints = {
        "kit_type": "molecular",
        "pathogen_type": "retrovirus",
        "pathogen_genome": "RNA",
        "identification_method": "sequencing_of_genetic_material",
        "detection_speed": "quick"
    }

    # --- Define Properties of Each Option from the Question ---
    # This maps the option letter to its described methodology.
    options_properties = {
        'A': {
            "kit_type": "serological",
            "identification": "IgG antibodies",
            "technique": "ELISA",
            "is_quick_for_early_detection": False,
            "reasoning": "This is a serological test, not molecular. IgG is a late-stage antibody, making it unsuitable for quick, early diagnosis."
        },
        'B': {
            "kit_type": "molecular",
            "identification": "cDNA sequencing",
            "technique": "real-time PCR",
            "is_quick_for_early_detection": True,
            "reasoning": "This is the correct workflow. A retrovirus (RNA genome) requires conversion to cDNA for sequencing. Real-time PCR is a quick and accurate molecular detection method."
        },
        'C': {
            "kit_type": "molecular",
            "identification": "symptoms",
            "technique": "nested PCR",
            "is_quick_for_early_detection": False, # Nested PCR is slower than real-time
            "reasoning": "Identifying a virus based on symptoms is scientifically invalid for designing a specific molecular test."
        },
        'D': {
            "kit_type": "molecular",
            "identification": "DNA sequencing",
            "technique": "PCR",
            "is_quick_for_early_detection": True,
            "reasoning": "This is incorrect because a retrovirus has an RNA genome. Direct DNA sequencing is the wrong first step and misses the crucial reverse transcription to cDNA."
        }
    }

    # --- Determine the correct option based on logic ---
    correct_option = None
    for option, props in options_properties.items():
        # Check 1: Must be a molecular kit
        if props["kit_type"] != constraints["kit_type"]:
            continue
        
        # Check 2: Identification method must be valid for an RNA virus
        # 'cDNA sequencing' is valid. 'DNA sequencing' and 'symptoms' are not.
        if props["identification"] == "DNA sequencing":
            continue
        if props["identification"] == "symptoms":
            continue

        # Check 3: Must be suitable for quick detection
        if not props["is_quick_for_early_detection"]:
            continue
            
        # If all checks pass, this is the correct option
        correct_option = option
        break

    # --- Extract the LLM's chosen answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>> in the provided text."
    
    llm_choice = match.group(1)

    # --- Compare LLM's choice with the logically derived correct answer ---
    if llm_choice == correct_option:
        return "Correct"
    else:
        incorrect_reason = options_properties[llm_choice]["reasoning"]
        correct_reason = options_properties[correct_option]["reasoning"]
        return (f"Incorrect. The chosen answer '{llm_choice}' is wrong. "
                f"Reason: {incorrect_reason}\n"
                f"The correct answer is '{correct_option}'. "
                f"Reason: {correct_reason}")


# In a real scenario, you would pass the full question and answer text.
# For this example, we'll use the provided final answer block.
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
There has been an outbreak of an viral infectious disease in your city. You have been asked to design a molecular diagnostic kit for quick detection of this retrovirus. How would you go about designing this kit?

A) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick 
 and accurate diagnosis.
B) First you identify the virus by performing cDNA sequencing.  You then develop a real time PCR kit for quick and accurate diagnosis.
C) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.
D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for  quick and accurate diagnosis
"""

llm_answer = """
Here is a step-by-step analysis of the question to determine the correct answer.

**1. Deconstruct the Question's Core Requirements:**
*   **Task:** Design a **molecular diagnostic kit**. This is a critical term. It means the kit must detect the pathogen's genetic material (nucleic acids like DNA or RNA), not the host's immune response (like antibodies) or other proteins.
*   **Pathogen:** The target is a **retrovirus**. This is the most important biological detail. Retroviruses are RNA viruses. They possess a unique enzyme called reverse transcriptase, which converts their RNA genome into a DNA copy (cDNA) as part of their replication cycle.
*   **Goal:** The kit must allow for **quick detection**. In an outbreak, this means it must be sensitive enough to detect early-stage infections and provide results rapidly.

**2. Analyze Each Option Against the Requirements:**

*   **A) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Careful Point 1 (Molecular vs. Serological):** This describes a serological (or immunological) test, not a molecular one. It detects antibodies produced by the host's immune system, not the virus's genetic material. This fails the first requirement.
    *   **Careful Point 2 (Timing):** IgG antibodies are markers of a later-stage or past infection. They take weeks to develop to detectable levels. A test for IgG is therefore unsuitable for "quick detection" of an active, early infection during an outbreak.

*   **B) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Careful Point 1 (Handling a Retrovirus):** This option correctly addresses the RNA nature of a retrovirus. To sequence its genome, the RNA must first be converted to complementary DNA (cDNA). This is the correct initial step.
    *   **Careful Point 2 (Diagnostic Technology):** Real-time PCR (specifically, RT-qPCR, which includes the reverse transcription step) is the gold standard for modern molecular diagnostics. It is extremely fast ("quick"), highly sensitive, and specific ("accurate"). It is the ideal technology for this scenario. This option describes the correct, state-of-the-art workflow.

*   **C) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Careful Point 1 (Identification Method):** This premise is scientifically invalid. You cannot determine the precise genetic sequence of a virus from clinical symptoms (e.g., fever, cough), as these are often non-specific. To design a PCR test, you need the exact genetic sequence to create specific primers.
    *   **Careful Point 2 (Technology Choice):** While nested PCR is a molecular technique, it is generally slower and more prone to contamination than real-time PCR, making it less suitable for a "quick" high-throughput diagnostic kit.

*   **D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Careful Point 1 (Biological Error):** This option misunderstands the fundamental biology of a retrovirus. The virus particle itself contains an RNA genome, not a DNA genome. Attempting to perform direct *DNA* sequencing on the virus would fail. This option misses the essential reverse transcription step (RNA -> cDNA).

**Conclusion:**

Option B is the only one that presents a scientifically correct, logical, and modern workflow. It correctly identifies the need to convert the retroviral RNA to cDNA for sequencing and then proposes using the most appropriate technology (real-time PCR) for rapid and accurate molecular diagnosis during an outbreak.

<<<B>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)