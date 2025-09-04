import re

def check_answer_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function encodes the scientific principles of molecular diagnostics for retroviruses
    and evaluates the chosen option against these principles.

    Args:
        question (str): The question text (for context, not directly used in logic).
        llm_answer_text (str): The full text of the LLM's response, including the final answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Extract the final answer from the LLM's response (e.g., <<<A>>>)
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find the final answer in the format <<<X>>> in the provided text."
    
    llm_choice = match.group(1)

    # Define the scientific constraints based on the question
    constraints = {
        "task": "molecular_diagnostic_kit",
        "pathogen": "retrovirus",
        "pathogen_genome": "RNA",
        "requirement": "quick_and_accurate"
    }

    # Define the properties of each option
    options = {
        'A': {
            "identification": "cDNA_sequencing",
            "kit_type": "real_time_PCR",
            "is_molecular": True,
            "handles_rna_genome": True, # cDNA sequencing implies RNA -> cDNA conversion
            "is_quick_for_outbreak": True, # Real-time PCR is the gold standard for speed
            "is_scientifically_valid": True
        },
        'B': {
            "identification": "IgG_antibodies",
            "kit_type": "ELISA",
            "is_molecular": False, # ELISA for antibodies is a serological/immunological test
            "handles_rna_genome": False, # Does not target the genome
            "is_quick_for_outbreak": False, # IgG indicates a late or past infection, not suitable for early/quick diagnosis
            "is_scientifically_valid": True # The method itself is valid, but not for this specific task
        },
        'C': {
            "identification": "symptoms",
            "kit_type": "nested_PCR",
            "is_molecular": True,
            "handles_rna_genome": False, # Identification method is invalid
            "is_quick_for_outbreak": False, # Identification method is impossible
            "is_scientifically_valid": False # Cannot design a molecular test from symptoms
        },
        'D': {
            "identification": "DNA_sequencing",
            "kit_type": "PCR",
            "is_molecular": True,
            "handles_rna_genome": False, # Incorrectly assumes a DNA genome for a retrovirus
            "is_quick_for_outbreak": True,
            "is_scientifically_valid": False # The first step is biologically incorrect
        }
    }

    # Determine the correct option based on the constraints
    correct_option = None
    for option_key, properties in options.items():
        if (properties["is_molecular"] and
            properties["handles_rna_genome"] and
            properties["is_quick_for_outbreak"] and
            properties["is_scientifically_valid"]):
            correct_option = option_key
            break

    # Check if the LLM's choice matches the correct option
    if llm_choice == correct_option:
        return "Correct"
    else:
        # Provide a reason why the LLM's choice is wrong
        chosen_properties = options[llm_choice]
        
        if not chosen_properties["is_molecular"]:
            return f"Incorrect. The chosen answer '{llm_choice}' describes a serological/immunological test (ELISA for antibodies), but the question asks for a 'molecular diagnostic kit', which must target the pathogen's nucleic acids (RNA/DNA)."
        
        if not chosen_properties["handles_rna_genome"]:
            if chosen_properties["identification"] == "DNA_sequencing":
                return f"Incorrect. The chosen answer '{llm_choice}' is wrong because a retrovirus has an RNA genome. Direct 'DNA sequencing' is not the correct first step; the RNA must first be converted to cDNA."
            else: # This covers the 'symptoms' case
                 return f"Incorrect. The chosen answer '{llm_choice}' is wrong because its identification method is invalid. It is scientifically impossible to design a specific molecular test based on patient symptoms."

        if not chosen_properties["is_quick_for_outbreak"]:
            return f"Incorrect. The chosen answer '{llm_choice}' describes a method not suitable for 'quick' diagnosis in an outbreak. Detecting IgG antibodies is for late-stage or past infections, not early detection."
            
        if not chosen_properties["is_scientifically_valid"]:
             return f"Incorrect. The chosen answer '{llm_choice}' is based on a scientifically invalid premise. For example, you cannot design a PCR kit from symptoms, or directly DNA sequence an RNA virus."

        return f"Incorrect. The chosen answer '{llm_choice}' is not the best fit. The correct answer is '{correct_option}' because it correctly identifies the need to convert the retroviral RNA to cDNA for sequencing and then uses the gold-standard real-time PCR for quick and accurate diagnosis."


# The user's provided question and answer
question = """
There has been an outbreak of an viral infectious disease in your city. You have been asked to design a molecular diagnostic kit for quick detection of this retrovirus. How would you go about designing this kit?

A) First you identify the virus by performing cDNA sequencing.  You then develop a real time PCR kit for quick and accurate diagnosis.
B) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick 
 and accurate diagnosis.
C) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.
D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for  quick and accurate diagnosis
"""

llm_answer = """
Here is a step-by-step analysis of the question and the provided options to determine the most accurate answer.

**1. Deconstruct the Question**

*   **Core Task:** Design a **molecular diagnostic kit**. This specifies that the method must detect the genetic material (nucleic acids like DNA or RNA) of the pathogen. This immediately makes methods that detect proteins or antibodies (like ELISA) less likely to be the primary answer.
*   **Pathogen Type:** The disease is caused by a **retrovirus**. This is a critical detail. Retroviruses are a class of RNA viruses. Their defining feature is the enzyme reverse transcriptase, which they use to convert their RNA genome into a DNA copy (cDNA) inside a host cell.
*   **Kit Requirements:** The kit must provide **quick and accurate diagnosis**. This is essential for managing an outbreak, where speed and reliability are paramount to identify infected individuals and prevent further spread.

**2. Analyze Each Option Step-by-Step**

*   **A) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Step 1: Identify by cDNA sequencing.** Since the virus is a retrovirus (an RNA virus), its genetic material cannot be sequenced directly by standard DNA sequencers. The RNA must first be converted into complementary DNA (cDNA) using reverse transcriptase. Sequencing this cDNA is the correct and necessary first step to determine the virus's unique genetic code. This step is scientifically sound.
    *   **Step 2: Develop a real-time PCR kit.** Once the viral sequence is known, specific primers and probes can be designed to target unique regions of the viral genome. Real-time PCR (specifically, Real-Time Reverse Transcription PCR or RT-qPCR) is the gold standard for molecular diagnostics. It is highly sensitive (detects tiny amounts of virus), specific (doesn't cross-react with other viruses), and fast (provides results in a few hours). This perfectly matches the "quick and accurate" requirement.
    *   **Conclusion:** This option describes the most logical, modern, and scientifically correct workflow.

*   **B) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Careful Point 1 (Molecular vs. Serological):** This describes a serological (or immunological) test, not a molecular one. It detects the host's immune response (antibodies), not the virus's genetic material.
    *   **Careful Point 2 (Timing):** IgG antibodies are part of the long-term immune response. They take several days to weeks to become detectable after infection. This "window period" means the test is not suitable for *quick* diagnosis of an active, early-stage infection, which is crucial during an outbreak.
    *   **Conclusion:** This option is incorrect because it's not a molecular test and it's not suitable for early detection.

*   **C) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Careful Point 1 (Identification Method):** Identifying a virus based on symptoms is scientifically impossible for the purpose of designing a molecular test. Symptoms like fever, cough, and fatigue are non-specific and can be caused by countless different pathogens. You cannot derive a genetic sequence from a list of symptoms.
    *   **Careful Point 2 (Test Design):** It is impossible to design the specific primers needed for any type of PCR (including nested PCR) without first knowing the target genetic sequence.
    *   **Conclusion:** The premise of this option is fundamentally flawed.

*   **D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Careful Point (Genetic Material):** This is incorrect because a retrovirus particle contains an RNA genome, not a DNA genome. Attempting to perform direct DNA sequencing on the isolated virus would fail. This option misses the essential reverse transcription step (converting RNA to cDNA).
    *   **Conclusion:** This option is incorrect because it misunderstands the basic biology of a retrovirus.

**3. Final Synthesis**

After a thorough analysis, Option A is the only one that presents a correct and logical workflow. It correctly accounts for the RNA nature of the retrovirus by including the cDNA step and proposes using the most appropriate modern technology (real-time PCR) for rapid and accurate molecular diagnosis during an outbreak. The other options are flawed due to incorrect biological assumptions, inappropriate test types, or scientifically impossible methods.

<<<A>>>
"""

# Run the check
result = check_answer_correctness(question, llm_answer)
print(result)