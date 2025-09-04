import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the molecular diagnostic kit question.

    The function verifies the answer based on key scientific principles:
    1.  A 'molecular' kit detects nucleic acids (DNA/RNA).
    2.  A 'retrovirus' has an RNA genome.
    3.  Detecting/sequencing RNA requires a reverse transcription step to create cDNA.
    4.  A 'quick and accurate' kit for an outbreak implies a modern, fast, and sensitive method like real-time PCR.
    5.  Serological tests (like ELISA for IgG) are not molecular and are not ideal for early/quick diagnosis.
    6.  Symptoms are not a valid basis for designing a specific molecular test.
    """

    # Step 1: Extract the final answer choice from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the standard format '<<<X>>>'."
    
    selected_option = match.group(1)

    # Step 2: Define the analysis for each option based on scientific facts.
    # The correct option is 'A'.
    correct_option = 'A'

    analysis = {
        'A': {
            "is_correct": True,
            "reason": "This is the correct approach. A retrovirus is an RNA virus, so its genome must be converted to cDNA for sequencing. Real-time PCR is the gold standard for quick and accurate molecular diagnosis."
        },
        'B': {
            "is_correct": False,
            "reason": "This is incorrect. Identifying a virus based on symptoms is not a valid method for designing a specific molecular test, as the genetic sequence is required to design PCR primers."
        },
        'C': {
            "is_correct": False,
            "reason": "This is incorrect. An ELISA kit for IgG is a serological test, not a molecular one. Furthermore, IgG antibodies appear late in an infection, making this method unsuitable for 'quick' early diagnosis."
        },
        'D': {
            "is_correct": False,
            "reason": "This is incorrect. A retrovirus has an RNA genome, not a DNA genome. Direct DNA sequencing is the wrong initial step and misses the essential reverse transcription step to create cDNA."
        }
    }

    # Step 3: Compare the LLM's selected option with the correct analysis.
    if selected_option == correct_option:
        return "Correct"
    else:
        # The LLM chose an incorrect option. Provide the reason why it's wrong.
        reason_for_error = analysis[selected_option]['reason']
        return f"The answer is incorrect. The selected option '{selected_option}' is wrong because: {reason_for_error}"

# The provided answer from the LLM
llm_response = """
Here is a step-by-step analysis of the question and the provided options.

**1. Deconstruct the Question**
*   **Core Task:** Design a **molecular diagnostic kit**. This specifies that the kit must detect the genetic material (nucleic acids like DNA or RNA) of the pathogen. This excludes methods that detect the host's immune response (like antibodies) or viral proteins.
*   **Pathogen Type:** The disease is caused by a **retrovirus**. This is a critical detail. Retroviruses are RNA viruses. Their defining characteristic is that they use an enzyme called reverse transcriptase to create a DNA copy (complementary DNA or cDNA) of their RNA genome.
*   **Kit Requirements:** The kit must enable **quick and accurate** diagnosis. This is crucial for managing an outbreak, where identifying infected individuals early is necessary to prevent further spread. This points towards a method that is both fast and highly sensitive/specific.

**2. Analyze Each Option Step-by-Step**

*   **A) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Step 1: Identify by cDNA sequencing.** This is the correct initial step. Since the virus is a retrovirus (an RNA virus), its genetic material must first be converted to complementary DNA (cDNA) using reverse transcriptase before it can be sequenced. This allows for the precise identification of the virus's genetic code.
    *   **Step 2: Develop a real-time PCR kit.** This is the ideal second step. Real-time PCR (specifically, Real-Time Reverse Transcription PCR or RT-qPCR) is the gold standard for molecular diagnostics. It is highly sensitive, specific, quantitative, and provides results very quickly (within hours), perfectly matching the "quick and accurate" requirement.
    *   **Conclusion:** This option describes a logical, scientifically correct, and state-of-the-art approach.

*   **B) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Step 1: Identify by symptoms.** This is scientifically invalid. Symptoms like fever or cough are non-specific and can be caused by many different pathogens. You cannot determine the unique genetic sequence of a virus, which is necessary to design PCR primers, simply by observing a patient's symptoms.
    *   **Conclusion:** This option is fundamentally flawed from the very first step.

*   **C) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Step 1 & 2: Identify IgG and develop an ELISA kit.** This describes a serological test, not a *molecular* one, as it detects the host's immune response (antibodies) rather than the virus's genetic material. Furthermore, IgG antibodies are part of the long-term immune response and take several weeks to appear. This creates a "window period" where an infected person is contagious but would test negative, making it unsuitable for *quick*, early diagnosis during an outbreak.
    *   **Conclusion:** This option fails on two counts: it is not a molecular test, and it is not suitable for early detection.

*   **D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Step 1: Identify by DNA sequencing.** This is incorrect. The virus particle itself contains an RNA genome. Attempting to directly perform DNA sequencing on the isolated virus would fail to identify its genetic material. This option misses the critical reverse transcription step required to create cDNA from the viral RNA.
    *   **Conclusion:** This option is flawed because it misunderstands the basic biology of a retrovirus.

**Final Synthesis**

Comparing the four options, Option A is the only one that presents a correct and logical workflow. It correctly acknowledges the RNA nature of the retrovirus by starting with cDNA synthesis and sequencing. It then proposes using the most appropriate modern technology (real-time PCR) for the rapid and accurate molecular diagnosis required during an outbreak.

<<<A>>>
"""

# Run the check and print the result
result = check_answer_correctness(llm_response)
print(result)