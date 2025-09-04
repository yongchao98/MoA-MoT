import re

def check_answer(question, llm_answer):
    """
    Checks the correctness of the LLM's answer to the biology question.

    The logic is based on the following biological principles:
    1.  A "lineage-specific promoter" is inactive in undifferentiated iPSCs. The red signal only appears after differentiation begins.
    2.  Apoptosis (green TUNEL signal) is a common and expected event in early embryos, especially for injected cells that fail to integrate.
    3.  The most significant experimental finding would be the event that links differentiation (red) and cell fate (green), which is colocalization.
    """

    # --- Define Biological Principles ---
    principles = {
        "promoter_activity": "A lineage-specific promoter is OFF in undifferentiated iPSCs. Red signal requires differentiation, which takes time.",
        "apoptosis_signal": "A green TUNEL signal indicates apoptosis. Apoptosis is expected in the developing embryo and for many injected cells. Therefore, a green signal is expected.",
        "mRaspberry_location": "mRaspberry is a cytoplasmic protein unless specifically engineered with targeting sequences. The promoter does not control subcellular location.",
        "key_finding": "The experiment aims to link iPSC fate (differentiation) with apoptosis. The most significant finding that directly addresses this is the colocalization of the red (differentiation) and green (apoptosis) signals."
    }

    # --- Extract the choice from the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Incorrect: The answer is not in the required format '<<<A>>>', '<<<B>>>', etc."
    
    choice = match.group(1)

    # --- Evaluate the choice against the principles ---
    if choice == 'A':
        # A) green signal colocalizes with the red signal
        # This aligns with the "key_finding" principle. It describes a plausible and significant
        # biological event where a cell that starts to differentiate is then eliminated by apoptosis.
        # This is the most informative outcome among the choices.
        return "Correct"
    
    elif choice == 'B':
        # B) cell line-specific red signals label different organelles
        return f"Incorrect: This choice violates a key principle. Reason: {principles['mRaspberry_location']}"

    elif choice == 'C':
        # C) cytoplasmic localization of the red signal
        # While technically true that the protein is cytoplasmic, this is a static property of the tool,
        # not a dynamic experimental finding or the "first thing you notice" about the biological process.
        # The question asks for an observation on cell fate, not protein biochemistry.
        return "Incorrect: This describes a static property of the reporter protein, not a significant biological finding about cell fate, which is the focus of the experiment."

    elif choice == 'D':
        # D) there is no green signal
        return f"Incorrect: This choice is biologically improbable. Reason: {principles['apoptosis_signal']}"

    else:
        return "Incorrect: The provided choice is not a valid option (A, B, C, or D)."

# The user-provided question and the LLM's answer to be checked
question = """
You want to create a mouse embryonic chimera with induced pluripotent stem cells from somatic cells of various tissues. You are interested in the iPSCs fate in the embryo and their co-localization with apoptotic events. The iPSC cells are labeled with a mRaspberry fused under a lineage-specific promoter for each cell type. You inject the dedifferentiated cells into the 32-cell stage blastocyst and track them over 48. Every 12 h, you prepare embryos with the TUNEL - FITC staining and observe under the confocal microscope.

What is the first thing you notice?

A) green signal colocalizes with the red signal
B) cell line-specific red signals label different organelles
C) cytoplasmic localization of the red signal
D) there is no green signal
"""

llm_answer = """
Here is a step-by-step analysis to determine the most likely answer.

### Step 1: Deconstruct the Experimental Setup

*   **Cells:** Induced pluripotent stem cells (iPSCs). By definition, these are in an undifferentiated, pluripotent state.
*   **Host:** A 32-cell stage mouse blastocyst.
*   **Red Label (for iPSCs):** The red fluorescent protein mRaspberry is used. Crucially, its gene is controlled by a **lineage-specific promoter**. This means the promoter is only activated, and the red protein is only produced, *after* the iPSC has started to differentiate into a specific cell type (e.g., a neuron, muscle cell, etc.). Therefore, at the moment of injection, the iPSCs are not red.
*   **Green Label (for Apoptosis):** A TUNEL-FITC stain is used. The TUNEL assay detects DNA fragmentation, a hallmark of late-stage apoptosis (programmed cell death). FITC is a green fluorophore. Any cell undergoing apoptosis, whether it's an injected iPSC or a host embryo cell, will be stained green.
*   **Question:** The question asks for the "first thing you notice," implying the earliest significant observation after the experiment begins (e.g., at the 12-hour time point).

### Step 2: Analyze the Expected Biological Events

1.  **Fate of Injected iPSCs:** Creating chimeras is an inefficient process. A large percentage of injected stem cells fail to properly integrate into the host embryo's inner cell mass. These cells are recognized as abnormal or misplaced and are rapidly eliminated by the embryo's quality control mechanisms.
2.  **Mechanism of Elimination:** The primary mechanism for removing these unwanted cells during development is apoptosis. Therefore, it is highly expected that many of the injected iPSCs will undergo apoptosis.
3.  **Timing of Gene Expression vs. Apoptosis:**
    *   **Apoptosis (Green Signal):** This can happen very quickly. Cells that fail to integrate can be targeted for elimination within hours. Apoptosis is also a normal part of embryogenesis, so a green signal is almost guaranteed to be present from the first time point.
    *   **Differentiation (Red Signal):** This is a slower process. The iPSC must receive signals from the embryonic environment, commit to a specific lineage, activate the promoter, and then transcribe and translate enough mRaspberry protein to become visibly red.
    *   **The Key Event:** A very common outcome is that an iPSC begins to differentiate (turning red) but does so in the wrong location or into an inappropriate cell type. The embryo then identifies this misplaced, newly-differentiated cell and eliminates it via apoptosis (turning it green).

### Step 3: Evaluate the Options

Let's evaluate the choices based on this biological understanding. The question asks for the "first thing you notice," which should be interpreted as the first *significant biological finding* that addresses the experimental goals.

*   **D) there is no green signal:** This is highly improbable. Apoptosis is a fundamental process in normal embryonic development and is the expected fate for many of the injected cells. A green signal will almost certainly be present.
*   **B) cell line-specific red signals label different organelles:** This is biologically incorrect. The promoter controls *if* a protein is made, not *where* it goes within the cell. Unless specifically engineered with targeting sequences (which is not mentioned), mRaspberry is a cytoplasmic protein.
*   **C) cytoplasmic localization of the red signal:** While this is a true statement about the mRaspberry protein, it is a static property of the tool, not a dynamic experimental result. It's a baseline technical detail, not the primary biological finding about cell fate.
*   **A) green signal colocalizes with the red signal:** This describes the key biological event at the heart of the experiment. It signifies that an iPSC that has begun to differentiate (turning red) is being eliminated by the embryo via apoptosis (turning green). Observing this colocalization is the first direct evidence that answers the research question about the fate of differentiating iPSCs and their link to apoptosis. While some cells may die before differentiating (appearing only green), the colocalization is the most informative and significant initial finding.

### Conclusion

The most direct and important observation that addresses the experimental question is the elimination of iPSCs that have begun to differentiate. This would be visualized under the microscope as a cell that is simultaneously red (due to the lineage-specific promoter activating upon differentiation) and green (due to the TUNEL stain marking it for apoptosis). Therefore, the colocalization of the two signals is the most likely "first thing you notice."

<<<A>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)