import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the iPSC chimera question.

    The function codifies the key biological principles of the experiment:
    1.  A lineage-specific promoter is inactive in undifferentiated iPSCs, so no red signal is expected initially.
    2.  Apoptosis (green signal) is a normal and expected event in embryogenesis and for rejected cells.
    3.  The colocalization of red (differentiated) and green (apoptotic) signals is a key, significant biological finding.
    4.  The subcellular localization of the reporter protein is a static property, not a dynamic finding.

    Args:
        llm_answer_text: The full text of the LLM's response, ending with the answer.

    Returns:
        A string indicating "Correct" or a reason for the answer being incorrect.
    """
    # --- Step 1: Define the ground truth based on biological principles ---

    # Principle 1: Promoter activity and red signal
    # A lineage-specific promoter is OFF in undifferentiated iPSCs.
    # Therefore, a red signal only appears AFTER differentiation begins.
    red_signal_requires_differentiation = True

    # Principle 2: Apoptosis and green signal
    # Apoptosis is normal in embryos and a common fate for injected cells.
    # Therefore, a green signal is expected.
    green_signal_is_expected = True

    # Principle 3: Evaluating the options
    option_correctness = {
        'A': {
            'is_correct': False,
            'reason': "Incorrect. Option A ('cell line-specific red signals label different organelles') is biologically wrong. A promoter controls gene expression timing, not the subcellular localization of the resulting protein. mRaspberry is cytoplasmic unless specifically engineered otherwise."
        },
        'B': {
            'is_correct': False,
            'reason': f"Incorrect. Option B ('there is no green signal') is biologically wrong. Apoptosis is expected in a developing embryo and is a common fate for injected cells that fail to integrate. Therefore, a green signal is expected (Principle: green_signal_is_expected is {green_signal_is_expected})."
        },
        'C': {
            'is_correct': True,
            'reason': "Correct. This option describes the first significant biological finding that addresses the experiment's core question. It represents the case where an iPSC begins to differentiate (turning red) but is then eliminated by the embryo via apoptosis (turning green). This colocalization is a key and expected event in chimera formation."
        },
        'D': {
            'is_correct': False,
            'reason': "Incorrect. Option D ('cytoplasmic localization of the red signal') is a less optimal answer. While true that the protein is cytoplasmic, this is a static property of the reporter tool, not a dynamic biological finding about cell fate. Option C describes a significant interaction between the injected cells and the host embryo, which is the central point of the experiment."
        }
    }
    correct_answer = 'C'

    # --- Step 2: Extract the final answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find the final answer in the standard '<<<X>>>' format."
    
    llm_choice = match.group(1)

    # --- Step 3: Compare the LLM's answer with the ground truth ---
    if llm_choice == correct_answer:
        return "Correct"
    else:
        return option_correctness.get(llm_choice, {}).get('reason', f"Incorrect. The provided answer '{llm_choice}' is not the best choice. The correct answer is '{correct_answer}'.")

# The final answer provided by the LLM to be checked
llm_response = """
Here is a step-by-step analysis to decide the final answer.

### Step 1: Deconstruct the Experimental Setup

1.  **Cells & Label:** Induced pluripotent stem cells (iPSCs) are used. They are undifferentiated. They are engineered to express a red fluorescent protein (mRaspberry).
2.  **Promoter - The Critical Detail:** The gene for the red protein is controlled by a **lineage-specific promoter**. This means the promoter is "off" when the cells are undifferentiated and only turns "on" when the cells begin to differentiate into a specific cell type.
3.  **Apoptosis Detection:** A TUNEL-FITC stain is used. This labels cells undergoing apoptosis (programmed cell death) with a green signal.
4.  **Procedure & Timeline:** The iPSCs are injected into an early-stage embryo (blastocyst). The system is observed every 12 hours. The question asks for the *first* thing you notice.

### Step 2: Analyze the Expected Biological Events

*   **Fate of the Red Signal:** Because the iPSCs are undifferentiated at injection and the red protein is controlled by a lineage-specific promoter, there should be **no red signal** at the very beginning of the experiment (e.g., at the 12-hour time point). The cells must first receive signals to differentiate, and this process takes time.
*   **Fate of the Green Signal:** Apoptosis (leading to a green signal) is a very common event in this type of experiment. It occurs for two main reasons:
    1.  It's a normal part of embryonic development to remove unneeded cells.
    2.  It's the primary way the host embryo eliminates injected cells that fail to integrate properly.
    Therefore, a green signal is **highly expected** from the start.
*   **The "Significant Finding":** The goal of the experiment is to link iPSC fate (differentiation) with apoptosis. A key event in chimera formation is the elimination of donor cells that differentiate inappropriately. The sequence would be:
    1.  iPSC starts to differentiate -> lineage-specific promoter activates -> cell turns **red**.
    2.  Embryo recognizes this differentiation as abnormal (wrong place/wrong time) -> triggers apoptosis -> cell turns **green**.
    This results in a cell that is both red and green.

### Step 3: Evaluate the Options and Synthesize the Answer

The question "What is the first thing you notice?" is slightly ambiguous.

*   **Literal Interpretation:** The very first observation in time (12h) would likely be some green cells (apoptosis) and no red cells (no differentiation yet). However, "no red signal" is not an option.
*   **"Key Finding" Interpretation:** The question is likely asking for the first *significant biological finding* that addresses the experimental goals. The death of undifferentiated cells (green only) is an expected artifact. The first observation that truly links differentiation fate with apoptosis is seeing a cell that is both red and green. This colocalization indicates that the process of differentiation itself is linked to the cells being eliminated.

Most of the provided candidate answers correctly identify this "key finding" interpretation as the most plausible one for a biology question. They reason that while some cells die before differentiating, the most important initial result is observing the fate of the cells that *do* begin to differentiate, which is often elimination.

Let's re-examine the options from the prompt:
*   A) cell line-specific red signals label different organelles: Incorrect. The promoter doesn't control protein location, and mRaspberry is cytoplasmic.
*   B) there is no green signal: Incorrect. Apoptosis is expected.
*   **C) green signal colocalizes with the red signal:** This matches the "key finding" interpretation. It describes the fate of iPSCs that attempt to differentiate but are then eliminated. This is a major, observable phenomenon in these experiments.
*   D) cytoplasmic localization of the red signal: While true that the signal would be cytoplasmic, this is a static property of the reporter protein, not a dynamic biological finding about cell fate.

**Conclusion:** Based on a careful analysis of the biological principles and the common interpretation of such experimental questions, the most significant and expected initial finding is the colocalization of the differentiation marker (red) and the apoptosis marker (green).

<<<C>>>
"""

# Run the check
result = check_answer_correctness(llm_response)
print(result)