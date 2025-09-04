import re

def check_correctness(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the retrovirus diagnostic kit question.

    This function programmatically evaluates the chosen option based on key scientific principles
    derived from the question:
    1. The kit must be 'molecular' (detecting nucleic acids).
    2. The pathogen is a 'retrovirus', which has an RNA genome. This means any DNA-based
       analysis (like sequencing or PCR) must first involve converting RNA to cDNA.
    3. The kit must be for 'quick' and 'accurate' diagnosis.
    """
    # Extract the single letter answer from the '<<<...>>>' format
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find the final answer in the required format '<<<...>>>'."
    
    chosen_option = match.group(1)

    # Define the properties and flaws of each option based on scientific facts
    options_analysis = {
        'A': {
            'is_correct': False,
            'reason': "The identification method is flawed. Designing a specific molecular test based on non-specific symptoms is not scientifically valid. The pathogen's genetic sequence is required."
        },
        'B': {
            'is_correct': True,
            'reason': "This option correctly follows the standard scientific procedure. For a retrovirus (RNA genome), the RNA is converted to cDNA for sequencing. A real-time PCR kit is the gold standard for quick, accurate, and sensitive molecular diagnosis."
        },
        'C': {
            'is_correct': False,
            'reason': "This fails on two counts. First, it describes an immunological test (ELISA), not a 'molecular' diagnostic kit. Second, detecting IgG antibodies is not suitable for 'quick' diagnosis of an early infection, as IgG appears late in the immune response."
        },
        'D': {
            'is_correct': False,
            'reason': "The identification method is incorrect. A retrovirus has an RNA genome, so direct 'DNA sequencing' would not work. This option misses the crucial and mandatory step of reverse transcription (converting RNA to cDNA)."
        }
    }

    # Check if the chosen option is the correct one
    if chosen_option in options_analysis:
        analysis_result = options_analysis[chosen_option]
        if analysis_result['is_correct']:
            return "Correct"
        else:
            # The chosen answer is incorrect. Return the reason why.
            return f"Incorrect. The chosen answer '{chosen_option}' is wrong. Reason: {analysis_result['reason']}"
    else:
        return f"Invalid option '{chosen_option}' was provided. Please choose from A, B, C, or D."

# The user-provided text containing the final answer to be checked.
llm_answer_to_check = """
Here is a step-by-step analysis of the question to determine the most accurate answer.

**Step 1: Deconstruct the Core Requirements of the Question**

1.  **Task:** Design a **molecular diagnostic kit**. This specifically refers to a test that detects the genetic material (nucleic acids like DNA or RNA) of the pathogen. This is a critical distinction from immunological tests that detect antibodies or antigens.
2.  **Pathogen:** The cause is a **retrovirus**. This is a key biological fact. Retroviruses are a type of RNA virus. Their genome is made of RNA, and they use an enzyme called reverse transcriptase to create a DNA copy (cDNA) of their genome as part of their replication cycle.
3.  **Kit Attributes:** The kit must be for **quick and accurate** detection. This is vital for managing an outbreak, implying the need for a method that is sensitive, specific, and has a fast turnaround time.

**Step 2: Evaluate Each Option Against the Core Requirements**

*   **A) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Careful Point:** Identifying a virus for the purpose of designing a molecular test based on symptoms alone is scientifically impossible. Different viruses can cause identical symptoms (e.g., fever, cough). To design a PCR test, you need the precise genetic sequence of the virus to create specific primers. This option's first step is fundamentally flawed.

*   **B) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Careful Point 1 (Identification):** Since the pathogen is a retrovirus (an RNA virus), its genetic material cannot be sequenced directly by standard DNA sequencers. The viral RNA must first be isolated and converted into complementary DNA (cDNA) using the reverse transcriptase enzyme. Sequencing this cDNA is the correct way to determine the virus's unique genetic code. This step is scientifically sound.
    *   **Careful Point 2 (Kit Development):** Real-time PCR (specifically, real-time RT-PCR, which incorporates the reverse transcription step) is the gold standard for modern molecular diagnostics. It is highly sensitive, specific, provides results rapidly ("quick"), and can be quantitative, making it extremely "accurate". This method is perfectly suited for an outbreak scenario. This option describes the correct, state-of-the-art workflow.

*   **C) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Careful Point 1 (Molecular vs. Immunological):** This describes an ELISA, which is an immunological test that detects the host's antibody response. It is not a *molecular* diagnostic test, failing a primary requirement of the question.
    *   **Careful Point 2 (Timing):** IgG antibodies are part of the long-term immune response and typically take weeks to become detectable after infection. This "window period" makes an IgG test unsuitable for *quick* or early diagnosis of an active infection, which is critical during an outbreak.

*   **D) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Careful Point:** This option is incorrect because it misunderstands the basic biology of a retrovirus. The virus particle contains an RNA genome, not a DNA genome. Attempting to perform direct DNA sequencing on the isolated virus would fail. This option misses the essential reverse transcription step.

**Conclusion**

After a thorough analysis, Option B is the only one that presents a scientifically valid and logical workflow that meets all the requirements of the question. It correctly identifies the need to convert the retroviral RNA to cDNA for sequencing and then proposes using the most appropriate modern technology (real-time PCR) for rapid and accurate molecular diagnosis.

<<<B>>>
"""

# Run the check
result = check_correctness(llm_answer_to_check)
print(result)