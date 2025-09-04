import re

def check_answer(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function evaluates the selected option based on key scientific principles:
    1.  A retrovirus has an RNA genome.
    2.  A "molecular" kit must target nucleic acids (RNA/DNA).
    3.  A "quick" kit for an outbreak requires a method suitable for early and rapid detection.
    """

    # Extract the single letter answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."
    
    selected_option = match.group(1)

    # Define the reasons why each option is correct or incorrect
    reasons = {
        'A': "Incorrect. Option A is wrong because identifying a virus based on symptoms is not specific or reliable enough to design a molecular test. A precise genetic sequence is required to create PCR primers.",
        'B': "Incorrect. Option B is wrong because a retrovirus has an RNA genome, not a DNA genome. The correct initial step is to convert the viral RNA to complementary DNA (cDNA) before sequencing, not to perform direct DNA sequencing.",
        'C': "Correct",
        'D': "Incorrect. Option D is wrong for two main reasons: 1) It describes an immunological test (ELISA), not a 'molecular' diagnostic kit as required by the question. 2) It targets IgG antibodies, which appear late in an infection and are therefore not suitable for 'quick' early diagnosis during an outbreak."
    }

    # The question text in the prompt has the options in a different order than the analysis.
    # Let's map the provided analysis text to the actual option letters.
    # A) symptoms -> nested PCR
    # B) DNA sequencing -> PCR
    # C) cDNA sequencing -> real time PCR
    # D) IgG antibodies -> ELISA
    
    # The logic in the `reasons` dictionary correctly maps to the options as defined above.
    # For example, the prompt's option A is about symptoms, and the code's reason for 'A' addresses that.
    # The prompt's option D is about IgG, and the code's reason for 'D' addresses that.
    
    return reasons.get(selected_option, f"Invalid option '{selected_option}' provided.")

# Example usage with the provided context:
# The user's provided analysis correctly identifies 'C' as the right answer.
# Let's assume the LLM's final output reflects this.

question_text = """
There has been an outbreak of an viral infectious disease in your city. You have been asked to design a molecular diagnostic kit for quick detection of this retrovirus. How would you go about designing this kit?

A) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.
B) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for  quick and accurate diagnosis
C) First you identify the virus by performing cDNA sequencing.  You then develop a real time PCR kit for quick and accurate diagnosis.
D) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick 
 and accurate diagnosis.
"""

llm_answer_text = """
Here is a step-by-step analysis of the question to determine the correct answer.

**1. Deconstruct the Core Task**

The question asks for the best method to design a **molecular diagnostic kit** for **quick detection** of a **retrovirus** during an outbreak.

*   **Careful Point 1: Retrovirus Biology.** The most critical piece of information is that the pathogen is a retrovirus. By definition, retroviruses have an RNA genome. To replicate, they use an enzyme called reverse transcriptase to convert their RNA into a DNA copy (cDNA). Any molecular technique aiming to sequence or amplify the virus's genetic code must account for this.

*   **Careful Point 2: "Molecular" Diagnostic Kit.** The request is specifically for a *molecular* kit. This means the kit must detect the virus's nucleic acids (RNA or DNA), not the host's immune response (like antibodies) or viral proteins.

*   **Careful Point 3: "Quick Detection".** In an outbreak, speed is essential for diagnosis to control the spread. The chosen method must be suitable for early detection (when viral load might be low but the person is infectious) and provide rapid results.

**2. Evaluate Each Option Based on the Careful Points**

*   **A) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Flaw:** The initial step is scientifically invalid. You cannot determine the precise genetic sequence of a virus, which is required to design PCR primers, based on clinical symptoms (e.g., fever, cough) as they are often non-specific.

*   **B) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Flaw:** This is incorrect because a retrovirus has an RNA genome, not a DNA genome. Attempting to perform direct DNA sequencing on the virus particle would fail. This option misses the essential reverse transcription step.

*   **C) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Correctness:** This option correctly follows the necessary scientific workflow.
        *   **Step 1 (cDNA sequencing):** It correctly acknowledges the RNA nature of the virus by first converting it to complementary DNA (cDNA) for sequencing. This is the proper way to identify the unique genetic code of a new RNA virus.
        *   **Step 2 (Real-time PCR kit):** It proposes developing a real-time PCR kit. This is the modern gold standard for rapid, highly sensitive, and specific molecular diagnostics, perfectly matching the "quick and accurate" requirement.

*   **D) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Flaw 1:** This describes an immunological (serological) test, not a *molecular* one, failing a primary requirement of the question.
    *   **Flaw 2:** IgG antibodies are produced late in the immune response, often weeks after infection starts. A test for IgG is therefore unsuitable for *quick, early* diagnosis during an outbreak, as it would miss currently infected individuals in the early stages.

**Conclusion**

Option C is the only choice that presents a scientifically sound, modern, and complete workflow that addresses all constraints of the question. It correctly handles the RNA nature of the retrovirus and proposes the most appropriate technology for rapid and accurate molecular diagnosis.

<<<C>>>
"""

# Run the check
result = check_answer(question_text, llm_answer_text)
print(result)