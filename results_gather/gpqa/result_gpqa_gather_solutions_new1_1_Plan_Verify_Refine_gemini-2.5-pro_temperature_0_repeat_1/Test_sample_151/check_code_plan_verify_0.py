import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer to a biology question about protein complexes.

    The function encodes the biological principles of the experiment described in the question
    to determine the correct answer and compares it with the provided answer.
    """

    # --- 1. Define the biological context and principles ---
    # The question describes an experiment on a yeast cell forming a "shmoo".
    # Key facts about this state:
    # - The cell is arrested in the G1 phase of the cell cycle.
    # - There is massive transcriptional activity to express mating-related genes.
    # - The cell is actively prevented from entering S phase, so DNA replication is inhibited.
    # The experiment (ChIP-MS) targets "active chromatin", which are regions of high transcription.

    # --- 2. Model the abundance of each complex based on the context ---
    # We assign a conceptual "activity/relevance" score. A lower score means less relevant
    # to the primary process (active transcription) and thus less abundant in the pulldown.
    complex_relevance = {
        "nucleosome histone complex": 3,  # Structural backbone of ALL chromatin. Ubiquitous and highly abundant.
        "enhancer protein complex": 2,    # Core machinery for the active process (transcription). Highly abundant.
        "pre-initiation complex": 2,      # Core machinery for the active process (transcription). Highly abundant.
        "pre-replication complex": 1      # Associated with an INHIBITED process (replication). Least abundant.
    }

    # --- 3. Determine the correct answer based on the model ---
    # The question asks for the LEAST observed complex.
    correct_complex_name = min(complex_relevance, key=complex_relevance.get)

    # Map the options from the question text to their names.
    options_mapping = {
        "A": "nucleosome histone complex",
        "B": "pre-replication complex",
        "C": "enhancer protein complex",
        "D": "pre-initiation complex"
    }

    correct_option_letter = None
    for letter, name in options_mapping.items():
        if name == correct_complex_name:
            correct_option_letter = letter
            break

    if correct_option_letter is None:
        return "Error in checker: Could not determine the correct option letter."

    # --- 4. Extract the LLM's final answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."
    
    llm_option_letter = match.group(1)

    # --- 5. Compare and generate the result ---
    if llm_option_letter == correct_option_letter:
        return "Correct"
    else:
        reason = (
            f"The provided answer is '{llm_option_letter}', which corresponds to the '{options_mapping.get(llm_option_letter)}'.\n"
            f"However, the correct answer is '{correct_option_letter}', which is the '{correct_complex_name}'.\n"
            f"Constraint Violated: The question asks for the LEAST observed complex in an assay for ACTIVE chromatin in a G1-arrested cell.\n"
            f"Reasoning: The cell is actively transcribing but NOT replicating. "
            f"The pre-replication complex is associated with the inhibited process of DNA replication. "
            f"All other complexes (nucleosome, enhancer, pre-initiation) are either structural components of all chromatin or are directly involved in the highly active process of transcription. "
            f"Therefore, the pre-replication complex would be the least observed."
        )
        return reason

# The user-provided question and the final answer from the LLM
question = """
You identified a new quorum-sensing peptide in the bacteria Lactobacillus acidophilus. You are highly interested if its role in intercellular communication is conserved among the eukaryotes. You treat yeast Saccharomyces cerevisiae with the peptide and interestingly observe the formation of the shmoo. You want to learn about the proteome of the active chromatin in the shmoo and recover it by chromatin immunoprecipitation followed by mass spectrometry. The least of proteins from which complex below will you observe in your assay?

A) nucleosome histone complex
B) pre-replication complex
C) enhancer protein complex
D) pre-initiation complex
"""

llm_answer_text = """
Here is a step-by-step analysis to decide the final answer.

**1. Deconstruct the Experimental Scenario and Biological State**

*   **Organism & Condition:** The experiment uses yeast (*Saccharomyces cerevisiae*) that has been induced to form a "shmoo."
*   **Key Biological Process:** Shmoo formation is a response to mating pheromones. A critical and defining feature of this response is that the yeast cell **arrests its cell cycle in the G1 phase**.
*   **Cellular Activity during G1 Arrest:** The cell is not preparing to divide. Instead, it undergoes a massive **transcriptional reprogramming**, activating hundreds of genes required for mating. Concurrently, the progression into S phase (where DNA replication occurs) is **actively inhibited**.
*   **Experimental Goal:** The experiment aims to recover the proteome (the set of proteins) from **active chromatin**. In this context, "active" primarily refers to regions undergoing high levels of gene transcription.

**2. Analyze Each Protein Complex in this Context**

*   **A) nucleosome histone complex:** Nucleosomes are the fundamental structural units of all chromatin, both active and inactive. They are the most basic and ubiquitous protein-DNA complexes. Any experiment that isolates chromatin will recover a very large quantity of histones. They will be highly abundant.

*   **C) enhancer protein complex & D) pre-initiation complex:** These two complexes are the core machinery of active gene transcription. Enhancer complexes bind to DNA to boost gene expression, and the pre-initiation complex (which includes RNA polymerase) assembles at gene promoters to start transcription. Since the shmoo response is defined by a massive wave of transcription, these complexes will be abundant on the active chromatin being studied.

*   **B) pre-replication complex (pre-RC):** This complex assembles on DNA at origins of replication during the G1 phase. Its function is to "license" the DNA, making it ready for replication in the S phase. While the shmoo is arrested in G1 (the correct phase for pre-RC assembly), the cell is actively prevented from entering S phase. Therefore, the function of the pre-RC is stalled. The cell's primary activity is transcription, not the inhibited process of replication.

**3. Compare Relative Abundance and Conclude**

The experiment is designed to find proteins on *active chromatin* in a cell that is actively transcribing but not replicating.

*   Complexes (C) and (D) are the engines of the cell's primary activity (transcription) and will be abundant.
*   Complex (A) is the fundamental structural component and will be the most abundant protein by mass.
*   Complex (B) is associated with a process (DNA replication) that is specifically inhibited in this cellular state. While present on the chromatin at specific sites (origins of replication), it is not involved in the "active" process of transcription that defines the shmoo state.

Therefore, the **pre-replication complex** is the least associated with the primary function of the chromatin in this state and would be the least observed complex in an assay targeting the proteome of active chromatin.

<<<B>>>
"""

# Run the check
result = check_answer(question, llm_answer_text)
print(result)