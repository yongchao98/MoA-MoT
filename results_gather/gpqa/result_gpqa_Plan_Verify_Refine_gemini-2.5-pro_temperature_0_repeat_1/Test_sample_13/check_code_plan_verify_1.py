import re

def check_correctness(llm_answer_text):
    '''
    This function checks the correctness of the provided answer to the ChIP-seq question.

    Question Analysis:
    - Experiment: ChIP-seq for IKAROS transcription factor in human B cells.
    - Condition 1: PFA fixation -> ChIP peaks are observed.
    - Condition 2: PFA+DSG fixation -> The same peaks disappear.
    - Question: Where are these disappearing peaks located?

    Biological Rationale for Checking:
    1.  PFA (Paraformaldehyde) is a short-range crosslinker known to cause artifacts by trapping proteins in highly crowded regions (even without stable binding).
    2.  PFA+DSG is a more stringent dual-crosslinking method that stabilizes protein complexes before linking to DNA, thereby reducing proximity-based artifacts.
    3.  The "disappearance" of peaks when using the more stringent method implies the original peaks were artifacts.
    4.  Such artifacts are known to occur at "hyper-ChIPable" regions, which are characterized by high protein density and open chromatin.
    5.  Active promoters and enhancers are the canonical examples of these regions. Therefore, option C is the correct answer.
    '''

    # The user's provided text is a reasoning process, not a final answer.
    # A correct final answer must be in the format <<<C>>>.
    # This code will check for the presence of such a formatted answer within the text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)

    if not match:
        # This handles the case where the provided text is just reasoning, like the example.
        return "The provided text is a reasoning process, not a final answer. It does not select one of the options A, B, C, or D in the required format '<<<X>>>'. Therefore, its correctness as a final answer cannot be checked."

    selected_option = match.group(1)

    if selected_option == 'C':
        return "Correct"
    elif selected_option == 'A':
        return "Incorrect. The disappearance of peaks with the more stringent PFA+DSG fixation points to an artifact common in regions of high protein density. While some transcription factors bind repeats, active promoters and enhancers are the canonical examples of regions prone to this specific artifact."
    elif selected_option == 'B':
        return "Incorrect. The signal is described as 'ChIP peaks', which are by definition non-random. The disappearance is a systematic effect, suggesting the peaks were artifacts at specific types of genomic locations, not random background noise."
    elif selected_option == 'D':
        return "Incorrect. While introns can contain enhancers, 'active promoters and enhancers' (Option C) is a more direct and accurate description of the regions known to be 'hyper-ChIPable' and prone to PFA-specific artifacts. The size of the gene is not the relevant factor."
    else:
        return "Invalid option selected."

# The user provided this text as the "answer" to be checked.
llm_response = '''
I understand the difference between PFA and PFA+DSG fixation. PFA is a short-range crosslinker for direct protein-DNA interactions, while PFA+DSG is a stronger, dual-crosslinking method that first links proteins into complexes (with DSG) and then links the complex to DNA (with PFA). The disappearance of peaks with the stronger PFA+DSG method suggests a specific type of experimental artifact. I will now research the function of the IKAROS transcription factor to understand where it typically binds and what kind of protein complexes it forms.

<search>IKAROS transcription factor function and binding sites in B cells</search>
'''

# Run the check
# print(check_correctness(llm_response))