import re

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer to the biological question.

    The logic follows these steps:
    1.  Analyze the core principle of the experiment: the use of a lineage-specific promoter.
    2.  Deduce the state of the red signal at the first time point (12h).
    3.  Evaluate each multiple-choice option based on this deduction.
    4.  Determine the correct answer by a process of elimination.
    5.  Compare the LLM's answer with the derived correct answer.
    """

    # --- Step 1 & 2: Analyze the experimental setup and deduce signal states ---
    # The most critical fact is that the red signal (mRaspberry) is controlled
    # by a lineage-specific promoter.
    # iPSCs are undifferentiated at the start.
    # Differentiation takes time, certainly more than 12 hours to produce a visible signal.
    # Therefore, at the first observation point (12h), there should be no red signal.
    red_signal_present_at_12h = False

    # --- Step 3: Evaluate the options based on the deduction ---
    # We create a dictionary to hold the reasons why each option is invalid.
    # If an option is valid, its value will be None.
    option_validity = {
        'A': None,
        'B': None,
        'C': None,
        'D': None,
    }
    
    # The question text has a different lettering for the options. Let's map them.
    # Question: A) there is no green signal, B) cell line-specific red signals..., C) green signal colocalizes..., D) cytoplasmic localization...
    # LLM Answer: A) there is no green signal, B) cell line-specific red signals..., C) green signal colocalizes..., D) cytoplasmic localization...
    # The mapping is direct.

    # Option B: "cell line-specific red signals label different organelles"
    if red_signal_present_at_12h is False:
        option_validity['B'] = "This option requires a red signal to be present, but the lineage-specific promoter would be inactive at 12h."
    
    # Option C: "green signal colocalizes with the red signal"
    if red_signal_present_at_12h is False:
        option_validity['C'] = "This option requires a red signal for colocalization, but the lineage-specific promoter would be inactive at 12h."

    # Option D: "cytoplasmic localization of the red signal"
    if red_signal_present_at_12h is False:
        option_validity['D'] = "This option requires a red signal to be localized, but the lineage-specific promoter would be inactive at 12h."

    # --- Step 4: Determine the correct answer by elimination ---
    possible_answers = [opt for opt, reason in option_validity.items() if reason is None]
    
    if len(possible_answers) != 1:
        # This case should not be reached with this logic, but it's good practice.
        return f"Logic Error: The analysis resulted in {len(possible_answers)} possible answers: {possible_answers}."
    
    correct_answer = possible_answers[0]

    # --- Step 5: Compare LLM's answer with the derived correct answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>>."

    llm_choice = match.group(1)

    if llm_choice == correct_answer:
        return "Correct"
    else:
        reason_for_llm_error = option_validity.get(llm_choice, "The chosen option is incorrect for other reasons.")
        return (f"Incorrect. The provided answer is <<<{llm_choice}>>>, but the correct answer is <<<{correct_answer}>>>. "
                f"The reasoning is as follows: The key constraint is that the red signal is controlled by a lineage-specific promoter. "
                f"In undifferentiated iPSCs at the 12-hour mark, this promoter would be inactive, meaning there is no red signal. "
                f"This fact invalidates options B, C, and D, as they all presuppose the existence of a red signal. "
                f"By elimination, A is the only possible answer. The chosen option <<<{llm_choice}>>> is wrong because: {reason_for_llm_error}")


# The user's provided question and the final answer from the LLM
question_text = """
You want to create a mouse embryonic chimera with induced pluripotent stem cells from somatic cells of various tissues. You are interested in the iPSCs fate in the embryo and their co-localization with apoptotic events. The iPSC cells are labeled with a mRaspberry fused under a lineage-specific promoter for each cell type. You inject the dedifferentiated cells into the 32-cell stage blastocyst and track them over 48. Every 12 h, you prepare embryos with the TUNEL - FITC staining and observe under the confocal microscope.

What is the first thing you notice?

A) there is no green signal
B) cell line-specific red signals label different organelles
C) green signal colocalizes with the red signal
D) cytoplasmic localization of the red signal
"""

llm_final_answer = """
Here is a step-by-step analysis of the problem to determine the final answer.

### Plan
1.  Deconstruct the experimental setup, focusing on the key components: the type of cells (iPSCs), the reporter system (mRaspberry under a lineage-specific promoter), the apoptosis marker (TUNEL), and the timeline ("first thing you notice" at 12h).
2.  Analyze the expected state of the red fluorescent signal at the first time point based on the principles of gene expression.
3.  Use the conclusion about the red signal to evaluate the feasibility of each multiple-choice option.
4.  Synthesize the findings to select the most logically sound answer by a process of elimination.

### Execution

**Step 1: Deconstruct the Experimental Setup**
*   **Cells:** Induced pluripotent stem cells (iPSCs) are used. By definition, these are in an undifferentiated state at the start of the experiment.
*   **Red Reporter:** The red fluorescent protein, mRaspberry, is controlled by a **lineage-specific promoter**. This is the most critical piece of information.
*   **Green Marker:** The TUNEL-FITC stain detects late-stage apoptosis (programmed cell death).
*   **Timeline:** The question asks for the "first thing you notice," with observations starting at 12 hours post-injection.

**Step 2: Analyze the Red Signal**
A lineage-specific promoter is a genetic switch that is only activated when a cell differentiates into a specific, specialized cell type (e.g., a muscle-specific promoter is active only in muscle cells). Since the iPSCs are injected in an undifferentiated state, this promoter will be inactive. The process of differentiation, followed by the transcription and translation of enough mRaspberry protein to be visible, takes a significant amount of time, almost certainly longer than 12 hours.
*   **Conclusion:** At the first observation point (12 hours), there should be **no red signal**.

**Step 3: Evaluate the Options Based on the "No Red Signal" Conclusion**
The fact that there is no red signal at the initial time point allows us to definitively evaluate three of the four options.

*   **B) cell line-specific red signals label different organelles:** This is impossible. There is no red signal to observe. (This option is also flawed because a promoter does not determine a protein's subcellular location).
*   **C) green signal colocalizes with the red signal:** This is impossible. Since there is no red signal, there can be no colocalization with it.
*   **D) cytoplasmic localization of the red signal:** This is impossible. There is no red signal to be localized.

**Step 4: Conclude by Elimination**
The fundamental principle of the experimental design—the inactivity of a lineage-specific promoter in an undifferentiated cell—logically invalidates options B, C, and D as the *first* observation. By a clear process of elimination, this leaves only option A.

*   **A) there is no green signal:** While apoptosis (green signal) is a common and expected event in early development and in response to injection, it is not a logical certainty. The TUNEL assay detects late-stage apoptosis, and it is plausible that at the very early 12-hour mark, widespread, detectable apoptosis has not yet occurred. More importantly, the other three options are rendered impossible by the core premise of the reporter system. Therefore, this must be the intended answer. The question is primarily designed to test the understanding that a lineage-specific promoter is silent in an undifferentiated cell.

<<<A>>>
"""

# Run the check
result = check_correctness(question_text, llm_final_answer)
print(result)