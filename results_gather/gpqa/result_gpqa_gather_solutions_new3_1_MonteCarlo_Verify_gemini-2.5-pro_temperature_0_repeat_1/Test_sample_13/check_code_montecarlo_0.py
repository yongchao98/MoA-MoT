import re

def check_correctness():
    """
    This function checks the correctness of the provided answer to a molecular biology question about ChIP-seq.

    Question: ChIP-seq peaks for the IKAROS transcription factor are seen with PFA fixation but disappear with PFA+DSG fixation. Where are these disappearing peaks most likely located?
    Options: A) At active promoters and enhancers, B) In the introns of large genes, C) At repeats, D) At random locations in the genome.
    Provided Answer: <<<C>>>

    The function evaluates two main competing hypotheses:
    1.  **Epitope Masking (favors Option A):** The stronger PFA+DSG fixation blocks the antibody binding site (epitope) at dense protein hubs like active promoters and enhancers.
    2.  **Insolubility/Artifact Removal (favors Option C):** The stronger fixation causes dense heterochromatic regions (composed of repeats, where IKAROS is known to bind) to become insoluble and get discarded during sample prep.

    The function concludes that the second hypothesis is more compelling due to the specific biology of IKAROS and the physical mechanism of insolubility, which better explains a complete disappearance of peaks.
    """
    
    # The final answer from the LLM being evaluated.
    llm_answer_text = "<<<C>>>"

    # Parse the answer from the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        # This case handles if the provided answer format is wrong.
        return f"Invalid answer format. Expected format like <<<A>>>, but got {llm_answer_text}"
    
    provided_answer = match.group(1)
    
    # The correct answer based on scientific reasoning.
    correct_answer = "C"
    
    if provided_answer == correct_answer:
        return "Correct"
    else:
        reasoning = (
            "The provided answer is incorrect. The most plausible explanation for the disappearing peaks combines the specific biology of the IKAROS protein with the physical properties of different chromatin regions. "
            "1. **IKAROS Biology:** IKAROS is well-documented to bind to pericentromeric heterochromatin, which is composed of highly repetitive DNA (repeats). "
            "2. **Fixation Effects:** These heterochromatic regions are extremely dense. The stronger PFA+DSG dual-fixation can over-crosslink these regions, creating large, insoluble aggregates. "
            "3. **Experimental Outcome:** These insoluble aggregates are typically pelleted and discarded during the initial chromatin preparation steps, meaning they are physically removed from the experiment before the antibody pulldown. This leads to a complete 'disappearance' of the signal from these regions. "
            "While epitope masking at promoters (Option A) is a possibility, the PFA+DSG method is generally used to *improve* signal at such functional sites. The catastrophic loss of signal is better explained by the insolubility of dense heterochromatic repeats. Therefore, the correct answer is C."
        )
        return f"Incorrect. The provided answer is '{provided_answer}', but the correct answer is '{correct_answer}'.\nReason: {reasoning}"

# Execute the check and print the result.
result = check_correctness()
print(result)