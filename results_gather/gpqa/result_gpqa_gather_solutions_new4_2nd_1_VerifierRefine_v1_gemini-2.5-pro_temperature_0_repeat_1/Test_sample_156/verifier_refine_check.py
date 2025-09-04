import re

def check_correctness_of_answer(question: str, answer: str) -> str:
    """
    Checks the correctness of the provided answer for the molecular diagnostic kit question.

    The function validates the answer based on key scientific principles:
    1.  **Molecular Kit**: The test must detect nucleic acids (RNA/DNA), not antibodies.
    2.  **Retrovirus**: The pathogen has an RNA genome, requiring reverse transcription (to cDNA) for DNA-based methods.
    3.  **Quick & Accurate**: The method must be fast and sensitive for early detection in an outbreak.

    Args:
        question: The question text.
        answer: The candidate answer text from the LLM, including the final <<<X>>> tag.

    Returns:
        A string indicating "Correct" or the reason for the error.
    """

    # --- Step 1: Define the ground truth based on scientific principles ---

    # A dictionary to hold the evaluation of each option.
    # The 'correct' option must satisfy all constraints.
    analysis = {
        'A': {
            'is_correct': True,
            'reason': "This is the correct workflow. 'cDNA sequencing' correctly handles the RNA genome of a retrovirus. 'Real-time PCR' is the gold standard for a quick, accurate, and sensitive molecular kit."
        },
        'B': {
            'is_correct': False,
            'reason': "This option is incorrect. The identification method, 'using...symptoms', is scientifically invalid for designing a specific molecular test. The genetic sequence cannot be determined from clinical symptoms."
        },
        'C': {
            'is_correct': False,
            'reason': "This option is incorrect. It misunderstands the biology of a retrovirus, which has an RNA genome. Direct 'DNA sequencing' would fail as it misses the essential reverse transcription step (RNA to cDNA)."
        },
        'D': {
            'is_correct': False,
            'reason': "This option is incorrect as it fails two key constraints. First, an ELISA kit is serological, not 'molecular'. Second, targeting IgG antibodies is not suitable for 'quick' early diagnosis during an outbreak."
        }
    }
    
    correct_option_id = 'A'

    # --- Step 2: Parse the LLM's chosen answer ---

    match = re.search(r'<<<([A-D])>>>\s*$', answer)
    if not match:
        return "Incorrect: The answer format is invalid. The final answer tag, e.g., '<<<A>>>', was not found at the end of the response."
    
    llm_choice = match.group(1)

    # --- Step 3: Compare the LLM's choice to the ground truth ---

    if llm_choice == correct_option_id:
        return "Correct"
    else:
        # The LLM chose an incorrect option. Provide the specific reason why its choice is wrong.
        reason_for_error = analysis[llm_choice]['reason']
        return f"Incorrect: The final answer <<< {llm_choice} >>> is wrong. {reason_for_error}"

# The provided question and answer from the prompt
question_text = """
There has been an outbreak of an viral infectious disease in your city. You have been asked to design a molecular diagnostic kit for quick detection of this retrovirus. How would you go about designing this kit?

A) First you identify the virus by performing cDNA sequencing.  You then develop a real time PCR kit for quick and accurate diagnosis.
B) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.
C) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for  quick and accurate diagnosis
D) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick 
 and accurate diagnosis.
"""

candidate_answer_text = """
Here is a step-by-step analysis of the question to determine the correct answer.

**1. Deconstruct the Question's Core Requirements**

*   **Task:** The primary task is to design a **molecular diagnostic kit**.
    *   **Careful Point:** This is a critical keyword. A "molecular" kit is designed to detect the genetic material (nucleic acids, i.e., DNA or RNA) of a pathogen. This requirement immediately excludes methods that detect the host's immune response (like antibodies) or viral proteins.

*   **Pathogen Type:** The pathogen is a **retrovirus**.
    *   **Careful Point:** This is the most important biological detail. Retroviruses are a specific class of RNA viruses. Their genome is made of RNA. They possess a unique enzyme called reverse transcriptase, which they use to create a DNA copy of their RNA genome inside a host cell. This DNA copy is called complementary DNA (cDNA). Any molecular technique aimed at the virus's genetic material must account for its RNA nature.

*   **Key Requirements:** The kit must be for **quick and accurate detection** during an **outbreak**.
    *   **Careful Point:** "Quick" means the test must provide results rapidly to allow for timely public health interventions like isolation and contact tracing. "Accurate" means the test must have high sensitivity (correctly identifying infected individuals, even with low viral loads) and high specificity (correctly identifying uninfected individuals). For an outbreak, early detection is paramount.

**2. Step-by-Step Evaluation of Each Option**

*   **A) First you identify the virus by performing cDNA sequencing. You then develop a real time PCR kit for quick and accurate diagnosis.**
    *   **Analysis:** This option describes a scientifically sound and modern workflow.
        *   **Step 1 (Identification):** It correctly addresses the RNA nature of a retrovirus by proposing **cDNA sequencing**. This implies the necessary first step of converting the viral RNA to complementary DNA (cDNA) using reverse transcriptase, which is the correct way to sequence an RNA virus.
        *   **Step 2 (Kit Development):** It proposes developing a **real-time PCR kit**. Real-time PCR (specifically, RT-qPCR, which includes the reverse transcription step) is the gold standard for modern molecular diagnostics. It is extremely fast, highly sensitive, and specific, making it the ideal choice for "quick and accurate" detection during an outbreak.
    *   **Conclusion:** This option meets all the requirements of the task.

*   **B) First you identify the virus using the information about the symptoms, provided by the affected patients. You then design a nested PCR kit based on the information for quick diagnosis.**
    *   **Analysis:** This option is fundamentally flawed.
        *   **Step 1 (Identification):** Identifying a virus based on symptoms is scientifically impossible for the purpose of designing a molecular test. Symptoms like fever or cough are non-specific and can be caused by many different pathogens. You cannot determine the unique genetic sequence of a virus from symptoms, which is essential for designing the specific primers needed for PCR.
    *   **Conclusion:** This method is invalid from the start.

*   **C) First you identify the virus by performing DNA sequencing. You then develop a PCR kit for quick and accurate diagnosis.**
    *   **Analysis:** This option contains a critical biological error.
        *   **Step 1 (Identification):** It proposes direct **DNA sequencing**. This is incorrect because a retrovirus has an RNA genome. Attempting to perform DNA sequencing on the virus particle would fail. This option misses the essential reverse transcription step (RNA -> cDNA).
    *   **Conclusion:** This option is based on a misunderstanding of retrovirus biology.

*   **D) First you identify the IgG antibodies for the infection. You then develop an ELISA kit targeting IgG antibodies for quick and accurate diagnosis.**
    *   **Analysis:** This option fails on two major criteria.
        *   **Kit Type:** An ELISA kit for antibodies is a serological (immunological) test, not a **molecular** one. It detects the host's immune response, not the virus's genetic material.
        *   **Speed of Detection:** It targets **IgG antibodies**, which are part of the long-term immune response. They take days to weeks to become detectable after infection. This creates a significant "window period" where an infected person would test negative, making it unsuitable for **quick**, early diagnosis during an outbreak.
    *   **Conclusion:** This option is not a molecular kit and is not suitable for early detection.

**Final Decision**

Based on the comprehensive analysis, Option A is the only one that presents a scientifically correct, logical, and modern workflow. It correctly accounts for the RNA nature of the retrovirus and proposes the most appropriate technologies for both initial identification (cDNA sequencing) and the final diagnostic kit (real-time PCR) to meet the requirements of quick and accurate detection during an outbreak.

<<<A>>>
"""

# Execute the check
result = check_correctness_of_answer(question_text, candidate_answer_text)
print(result)