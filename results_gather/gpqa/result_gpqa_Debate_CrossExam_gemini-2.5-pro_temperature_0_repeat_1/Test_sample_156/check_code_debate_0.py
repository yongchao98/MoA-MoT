import re

def check_answer_correctness(llm_response: str):
    """
    Checks the correctness of the LLM's answer for the retrovirus diagnostic kit question.

    The function validates the chosen option based on key scientific principles:
    1. Retroviruses are RNA viruses, requiring reverse transcription for DNA-based analysis.
    2. The request is for a 'molecular' kit (detecting nucleic acids), not immunological.
    3. 'Quick' detection in an outbreak favors methods like real-time PCR.
    4. Identification of a new virus requires sequencing, not just symptom analysis.
    """
    try:
        # Extract the final answer, which is expected in the format <<<X>>>
        answer_match = re.search(r'<<<([A-D])>>>', llm_response)
        if not answer_match:
            return "The answer format is incorrect. It should end with '<<<X>>>' where X is A, B, C, or D."
        
        provided_answer = answer_match.group(1)
    except Exception as e:
        return f"An error occurred while parsing the answer: {e}"

    # Define the analysis for each option based on scientific facts
    analysis = {
        "A": {
            "is_correct": False,
            "reason": "Option A is incorrect because identifying a virus based on symptoms is not specific enough for designing a molecular test. Additionally, nested PCR is generally slower and more prone to contamination than real-time PCR, making it less ideal for 'quick' diagnosis in an outbreak."
        },
        "B": {
            "is_correct": False,
            "reason": "Option B is incorrect because an ELISA for IgG antibodies is an immunological test, not a 'molecular' one as required by the question. Furthermore, IgG antibodies appear late in an infection, so this method is not suitable for 'quick' early detection."
        },
        "C": {
            "is_correct": False,
            "reason": "Option C is incorrect because a retrovirus has an RNA genome. It is not possible to perform 'DNA sequencing' directly on the virus. The viral RNA must first be converted to complementary DNA (cDNA) through reverse transcription."
        },
        "D": {
            "is_correct": True,
            "reason": "This option correctly outlines the modern standard procedure. For an RNA virus (retrovirus), its genome must be converted to cDNA for sequencing. Sequencing is essential to identify a new virus and design specific primers. A real-time PCR kit is the gold standard for quick, sensitive, and accurate molecular diagnosis during an outbreak."
        }
    }

    # Check if the provided answer is the correct one
    if provided_answer == "D":
        # The LLM's answer is correct.
        return "Correct"
    else:
        # The LLM's answer is incorrect, provide the reason why.
        return analysis[provided_answer]["reason"]

# The response provided by the other LLM
llm_response_from_prompt = """
Proposer: The initial analysis points strongly to option D. To confirm this, I will conduct a search to verify the standard procedures for developing a diagnostic kit for a new retrovirus. The search should confirm that sequencing is necessary to identify the virus and that real-time PCR is the preferred method for rapid and accurate detection.

Skeptic: While option D seems plausible, we need to be certain about the specifics. Is "cDNA sequencing" the standard term for sequencing an RNA virus? Is "real-time PCR" definitively the best choice over other PCR methods like nested PCR for this specific scenario (quick detection in an outbreak)? We should also confirm that immunological methods like ELISA (Option B) are indeed unsuitable for *quick, early* detection. The search queries should address these points to eliminate any ambiguity between the options.

Decided action: The search queries have been formulated to address the points raised by both the Proposer and the Skeptic. The goal is to confirm the entire workflow from identification to the final kit design for a novel retrovirus during an outbreak. The search results will be used to definitively select the best option among A, B, C, and D.

**SELF-DEBATE**

**Proposer:** The search results confirm that for an RNA virus, like a retrovirus, the genetic material must first be converted from RNA to complementary DNA (cDNA) before it can be sequenced or amplified by PCR. This immediately invalidates option C, which suggests direct DNA sequencing of a retrovirus. The results also emphasize that for rapid and accurate detection during an outbreak, real-time RT-PCR (Reverse Transcription-Polymerase Chain Reaction) is the gold standard. It is faster, more sensitive, and less prone to contamination than conventional PCR methods. This supports the "real time PCR kit" part of option D.

**Skeptic:** I agree that option C is incorrect. However, what about the other options? Option A suggests identifying the virus based on symptoms and then using nested PCR. Is symptom-based identification ever reliable enough? And is nested PCR a viable alternative?

**Proposer:** Identifying a virus based on symptoms alone is not a reliable method for designing a specific molecular test. Different viruses can cause similar symptoms. For a molecular kit, you need precise genetic information. While nested PCR is very sensitive, it is slower and more susceptible to contamination, making it less suitable for a "quick" diagnosis in a large-scale outbreak compared to real-time PCR. Real-time PCR provides results in hours and is designed for high-throughput screening.

**Skeptic:** That addresses Option A. Now, let's consider Option B: identifying IgG antibodies and developing an ELISA kit. The question asks for a *molecular* diagnostic kit. Is an ELISA kit considered molecular? And is it suitable for *quick* detection?

**Proposer:** An ELISA for IgG antibodies is an immunological assay, not a molecular one. Molecular diagnostics target the pathogen's nucleic acids (DNA or RNA). Furthermore, IgG antibodies are part of the body's later immune response, appearing well after the initial infection. Therefore, an IgG-based ELISA is not suitable for quick, early detection of an active infection during an outbreak. It is more useful for determining if someone has had a past infection or has developed immunity. Direct detection of the viral RNA is necessary for early diagnosis.

**Skeptic:** This confirms that options A, B, and C have significant flaws. Option D proposes identifying the virus by performing cDNA sequencing and then developing a real-time PCR kit. The search results strongly support this pathway. To identify a new RNA virus, sequencing is essential to determine its unique genetic code. This sequence information is then used to design specific primers and probes for a real-time RT-PCR assay. This method is highly sensitive, specific, and provides rapid results, which is critical for managing an outbreak.

**Proposer:** Exactly. The entire process described in option D aligns perfectly with the modern, standard approach for developing a molecular diagnostic for a novel RNA virus outbreak. You first sequence the viral genome (via cDNA) to know what you're looking for, and then you develop a real-time PCR test to detect it quickly and accurately in patients.

**Convergence:** Both Proposer and Skeptic agree that Option D presents the most accurate and scientifically sound method. It correctly identifies the genetic material of a retrovirus (RNA), the necessary step to analyze it (cDNA synthesis), the method for identification (sequencing), and the best-suited technology for a rapid diagnostic kit (real-time PCR).

<<<D>>>
"""

# Run the check
result = check_answer_correctness(llm_response_from_prompt)
print(result)