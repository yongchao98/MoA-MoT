import re

def check_correctness_of_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function validates the answer based on key scientific principles:
    1. A retrovirus is an RNA virus.
    2. A "molecular" kit must detect nucleic acids (RNA/DNA).
    3. A "quick" diagnostic for an outbreak should detect active, early-stage infection.
    """
    
    # Extract the final letter choice from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>'."
    
    provided_answer = match.group(1)

    # Define the properties of each option based on the question
    options_analysis = {
        'A': {
            "is_correct": False,
            "reason": "It incorrectly suggests 'DNA sequencing' as the first step for a retrovirus, which is an RNA virus. The viral RNA must first be converted to cDNA."
        },
        'B': {
            "is_correct": False,
            "reason": "It describes an immunological test (ELISA for IgG), not a 'molecular' diagnostic kit. Furthermore, IgG antibodies are markers for late or past infection, not for 'quick' early diagnosis."
        },
        'C': {
            "is_correct": False,
            "reason": "It suggests identifying a virus based on 'symptoms' to design a PCR kit, which is scientifically impossible. A precise genetic sequence is required."
        },
        'D': {
            "is_correct": True,
            "reason": "This option correctly describes the standard modern workflow: 1) Convert viral RNA to cDNA for sequencing to identify the virus. 2) Develop a real-time PCR kit, which is the gold standard for quick and accurate molecular diagnosis."
        }
    }

    # Determine the correct option letter
    correct_option = None
    for option, details in options_analysis.items():
        if details["is_correct"]:
            correct_option = option
            break
            
    if correct_option is None:
        # This is an internal check, should not be triggered with current logic
        return "Error in checker: No correct option was defined."

    # Compare the provided answer with the correct answer
    if provided_answer == correct_option:
        return "Correct"
    else:
        reason_for_error = options_analysis[provided_answer]["reason"]
        return f"Incorrect: The provided answer is '{provided_answer}', but the correct answer is '{correct_option}'. The chosen answer is wrong because {reason_for_error}"

# The user's provided response to check
llm_response = """
**Proposer**: The most logical solution is outlined in option D. The core of the problem lies in identifying a **retrovirus** (an RNA virus) and creating a **quick, molecular** diagnostic kit. Option D correctly proposes converting the viral RNA to cDNA for sequencing, which is necessary to identify the virus's genetic code. It then suggests developing a real-time PCR kit, which is the modern standard for rapid, sensitive, and specific molecular detection.

**Skeptic**: I agree that D seems strong, but let's not dismiss the other options too quickly.
*   Option A suggests DNA sequencing. Could the virus's DNA form, once integrated into the host cell, be the target?
*   Option B suggests an ELISA for IgG. While maybe not for "early" detection, ELISA tests can be very fast to run. Is it completely unsuitable? The question asks for "quick detection," which could mean a fast test turnaround, not necessarily early-stage infection detection.
*   Option C mentions nested PCR. This is a highly sensitive molecular technique. Why is real-time PCR definitively better?
*   The core assumption is that we are dealing with a completely novel virus. If it's a known retrovirus family, couldn't we skip the full sequencing step?

**Proposer**: Let's address those points.
*   Targeting integrated proviral DNA (Option A's logic) is possible but less ideal for a *quick* outbreak diagnosis. The primary target for early detection is the circulating virus particles, which contain RNA. Detecting viral RNA in a sample like a nasal swab or blood is a more direct and typically more sensitive measure of active infection than searching for integrated DNA in host cells. Furthermore, the option states "performing DNA sequencing" to *identify* the virus, which is incorrect for the virus particle itself.
*   Regarding Option B, the question asks for a *molecular* kit. An ELISA is an immunological assay, not a molecular one. It detects proteins (antibodies), not nucleic acids. Furthermore, IgG antibodies are a marker of a later-stage or past infection, making them unsuitable for identifying newly infected individuals to control an outbreak.
*   Regarding Option C, while nested PCR is sensitive, it is more complex, has a higher risk of cross-contamination (leading to false positives), and is generally slower than real-time PCR, especially in a high-throughput setting needed for an outbreak. Real-time PCR is a "closed-tube" system, which minimizes contamination and provides results much faster, fitting the "quick" requirement better.
*   Even if the virus belongs to a known family, a new outbreak implies it could be a new strain with genetic variations. Sequencing is crucial to confirm its identity and ensure the primers and probes for the PCR test are specific to the circulating strain, guaranteeing accuracy.

**Skeptic**: Your points are valid. The distinction between a "molecular" kit and an "immunological" one is critical, ruling out Option B. The impracticality of using symptoms to design a PCR test rules out C. The biological inaccuracy of directly DNA sequencing an RNA virus rules out A. This leaves D as the only option that follows a scientifically sound, modern workflow for this specific task. The process of converting RNA to cDNA for sequencing and then using that data to design a real-time PCR test is the gold standard.

**Decided Action**: The debate has converged. The analysis of the core requirements—molecular, quick, accurate, for a retrovirus—systematically eliminates options A, B, and C. Option D describes the precise, standard, and most effective methodology used in modern molecular diagnostics for an RNA virus outbreak.

***

### Final Analysis

Here is a step-by-step breakdown of the reasoning to select the final answer.

1.  **Analyze the Core Task**: The request is to design a **molecular diagnostic kit**.
    *   **Careful Point**: This specifies that the test must detect the pathogen's genetic material (nucleic acids, i.e., RNA or DNA). This immediately makes Option B (developing an ELISA kit for IgG antibodies) incorrect, as that is an immunological test that detects the host's protein-based immune response, not the virus's nucleic acids.

2.  **Analyze the Pathogen**: The pathogen is a **retrovirus**.
    *   **Careful Point**: Retroviruses are RNA viruses. Their genetic material is RNA. Standard technologies for sequencing and amplification, like PCR, work on DNA. Therefore, a critical step is required: converting the viral RNA into complementary DNA (cDNA) using an enzyme called reverse transcriptase.
    *   This fact makes Option A ("identify the virus by performing DNA sequencing") incorrect. You cannot directly use DNA sequencing on the RNA genome of the virus particle.

3.  **Analyze the Identification Method**: A new virus in an outbreak must be identified before a specific test can be designed.
    *   **Careful Point**: Option C suggests identifying the virus based on patient symptoms. This is scientifically invalid for designing a molecular test. Symptoms like fever or cough are non-specific and could be caused by hundreds of different pathogens. To design a PCR test, you need the precise genetic sequence of the target virus to create specific primers, which cannot be determined from symptoms.

4.  **Analyze the Diagnostic Technology**: The kit must be for **quick and accurate** detection.
    *   **Careful Point**: Option D suggests developing a **real-time PCR kit**. Real-time PCR (specifically, real-time RT-PCR for an RNA target) is the modern gold standard for this purpose. It is extremely sensitive and specific ("accurate"), and it provides results in a few hours ("quick"). It is also less prone to contamination and more suitable for high-throughput screening during an outbreak than older methods like the nested PCR mentioned in Option C.

**Conclusion**:

By carefully considering each component of the question, Option D emerges as the only correct choice. It accurately describes the modern, standard workflow for developing a molecular diagnostic for a novel RNA virus:

1.  **cDNA sequencing**: Correctly handles the RNA nature of the retrovirus to determine its genetic sequence.
2.  **Real-time PCR kit**: Employs the best available technology for rapid, accurate, and scalable diagnosis during an outbreak.

<<<D>>>
"""

# Run the check
result = check_correctness_of_answer(llm_response)
print(result)