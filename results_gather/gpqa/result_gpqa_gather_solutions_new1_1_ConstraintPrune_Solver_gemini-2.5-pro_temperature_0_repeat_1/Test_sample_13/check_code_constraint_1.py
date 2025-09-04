import re

def check_chip_seq_answer(llm_answer_option: str):
    """
    Checks the correctness of an answer to the IKAROS ChIP-seq question.

    The question asks to identify the location of IKAROS ChIP-seq peaks that are
    present with PFA fixation but disappear with PFA+DSG fixation.

    The correct answer is B, "At repeats". This is based on two main lines of reasoning:
    1.  **Insolubility of Dense Chromatin:** IKAROS is known to bind pericentromeric
        heterochromatin, which is dense and repetitive. The strong PFA+DSG dual-crosslinking
        can make these regions insoluble and resistant to sonication, causing them to be
        lost during the ChIP procedure, thus making the peaks "disappear".
    2.  **Artifact Removal:** Repetitive regions are a known source of artifacts in ChIP-seq.
        PFA-only fixation might capture weak or non-specific interactions at these "sticky"
        regions. The more stringent PFA+DSG method either prevents this or boosts the signal
        from true binding sites so much that these weaker, artifactual peaks fall below the
        statistical detection threshold.

    Args:
        llm_answer_option (str): The letter option (A, B, C, or D) provided by the LLM.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    correct_option = 'B'
    
    # Normalize the input to handle formats like "<<<B>>>" or "B"
    match = re.search(r'[A-D]', llm_answer_option.upper())
    if not match:
        return f"Invalid answer format. Expected one of A, B, C, or D, but got '{llm_answer_option}'."
        
    submitted_option = match.group(0)

    if submitted_option == correct_option:
        return "Correct"
    
    error_reasons = {
        'A': "Incorrect. This option suggests the peaks are at active promoters and enhancers. While epitope masking at these sites is a possibility, a more specific and compelling explanation involves the known biology of IKAROS binding to repetitive DNA. Strong cross-linking of dense heterochromatic regions (repeats) can lead to insolubility and loss of signal, making 'At repeats' (B) the better answer.",
        'C': "Incorrect. ChIP-seq peaks, even artifactual ones, are typically not at truly random locations. The systematic disappearance of peaks points to a specific biochemical cause related to chromatin structure, not a random process.",
        'D': "Incorrect. 'In the introns of large genes' is too general. While introns can contain repeats or regulatory elements, the key feature causing the phenomenon is the nature of the chromatin (e.g., repetitive and dense), not simply its location within an intron. Option B is more specific and mechanistically relevant."
    }
    
    return error_reasons.get(submitted_option, "Invalid option provided.")

# Example usage:
# Let's check the various answers provided by the candidate LLMs.
# Many answers were confused about the option letters. We will check based on the letter they returned.
print(f"Checking Answer 3 (<<<B>>>): {check_chip_seq_answer('B')}")
print(f"Checking Answer 13 (<<<A>>>): {check_chip_seq_answer('A')}")
print(f"Checking Answer 1 (<<<C>>>): {check_chip_seq_answer('C')}")
print(f"Checking Answer 4 (<<<D>>>): {check_chip_seq_answer('D')}")