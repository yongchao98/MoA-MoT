import re

def check_mitochondria_experiment_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to a question about mitochondrial experiments.

    The function uses a "ground truth" dictionary based on established cell biology principles
    to determine which experiment is not helpful for studying mitochondrial function. It then
    compares this correct answer to the one provided by the LLM.

    Args:
        llm_answer_text: The full text of the LLM's response, which should include
                         a final answer in the format <<<X>>>.

    Returns:
        "Correct" if the LLM's answer is correct.
        A string explaining the error if the LLM's answer is incorrect or malformed.
    """
    # Ground truth based on established cell biology principles.
    # The question asks which experiment will NOT help.
    ground_truth = {
        'A': {
            'description': "Luciferase assay for ATP",
            'is_helpful': True,
            'reasoning': "This experiment measures ATP, the primary energy product of mitochondria. It is a valid, though sometimes indirect, way to assess mitochondrial function and its impact on cell viability."
        },
        'B': {
            'description': "Confocal microscopy with Mito-RTP stain",
            'is_helpful': True,
            'reasoning': "This experiment uses a mitochondria-targeted probe (Mito-RTP) to directly measure parameters like temperature or ROS, which reflect metabolic activity within the organelle."
        },
        'C': {
            'description': "Glucose uptake assay on isolated mitochondria",
            'is_helpful': False,
            'reasoning': "This experiment is fundamentally flawed. Mitochondria do not directly take up glucose; they import pyruvate, the product of glycolysis which occurs in the cytoplasm. Therefore, a glucose uptake assay on isolated mitochondria is biologically nonsensical."
        },
        'D': {
            'description': "Flow cytometry with JC-1 dye",
            'is_helpful': True,
            'reasoning': "This experiment uses the JC-1 dye to measure mitochondrial membrane potential, a critical and direct indicator of mitochondrial health and function."
        }
    }

    # Determine the correct answer from the ground truth.
    # The correct answer is the option that is NOT helpful.
    correct_option = None
    for option, details in ground_truth.items():
        if not details['is_helpful']:
            correct_option = option
            break

    # Extract the LLM's chosen answer from the text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>'."

    llm_choice = match.group(1)

    # Compare the LLM's choice with the correct answer and provide feedback.
    if llm_choice == correct_option:
        return "Correct"
    else:
        reason_for_correct_choice = ground_truth[correct_option]['reasoning']
        reason_for_llm_choice = ground_truth[llm_choice]['reasoning']
        
        return (f"Incorrect. The final answer should be {correct_option}, but the provided answer was {llm_choice}. "
                f"The question asks for the experiment that will NOT help. "
                f"The correct answer is {correct_option} because: {reason_for_correct_choice} "
                f"The chosen option {llm_choice} is incorrect because it describes a valid experiment: {reason_for_llm_choice}")

# The user-provided answer text to be checked
llm_answer = """
Here is a step-by-step analysis of the experimental options to determine which one will not help discover the mitochondrial role of the drug. The goal is to identify the experiment that is based on a flawed or irrelevant premise for studying mitochondria.

### Step-by-step Analysis of the Options

**A) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
*   **Principle:** This is a bioluminescence assay that measures the concentration of ATP. The enzyme luciferase uses ATP as a substrate to produce light, and the amount of light is proportional to the amount of ATP.
*   **Relevance:** Mitochondria are the primary producers of cellular ATP through oxidative phosphorylation. Therefore, measuring ATP levels is a fundamental way to assess the primary output of mitochondrial function. A drug that impairs mitochondria would almost certainly affect cellular ATP levels.
*   **Careful Point:** The experiment specifies measuring ATP in the "supernatant," which means it measures extracellular ATP. While the main ATP pool is intracellular, a significant release of ATP can be a strong indicator of cellular stress or cell death (necrosis), which could be triggered by a mitochondrially toxic drug. Therefore, this experiment is an indirect but still helpful method for screening for drug effects related to cellular energy and viability.

**B) Confocal fluorescence microscopy after Mito-RTP staining of the cells**
*   **Principle:** This experiment uses a fluorescent probe called "Mito-RTP." The "Mito-" prefix indicates that the probe is specifically designed to accumulate within the mitochondria. Probes of this type are used to measure specific parameters of mitochondrial function or status.
*   **Relevance:** Depending on the specific probe, Mito-RTP could measure mitochondrial temperature (a proxy for metabolic rate), reactive oxygen species (ROS), or membrane potential. All of these are direct indicators of mitochondrial activity and health. Confocal microscopy provides high-resolution images to visualize and quantify these changes specifically within the organelle.
*   **Careful Point:** This is a direct and targeted method to observe if the drug affects mitochondrial activity or induces stress. It is a valid and helpful experiment.

**C) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
*   **Principle:** This experiment involves two steps: first, isolating mitochondria from the rest of the cell, and second, performing a glucose uptake assay on this isolated mitochondrial fraction.
*   **Relevance:** This experiment is designed to measure the rate of glucose transport into the isolated mitochondria.
*   **Careful Point (The Fundamental Flaw):** This experimental design is based on a false biological premise. The process of glucose uptake from the environment is carried out by glucose transporters (like GLUTs) located on the **cell's outer plasma membrane**, not on the mitochondrial membranes. The initial breakdown of glucose (glycolysis) occurs in the **cytoplasm**, converting glucose into pyruvate. It is **pyruvate**, not glucose, that is then transported into the mitochondria to fuel the Krebs cycle. Isolated mitochondria lack the necessary machinery to take up glucose. Therefore, applying a glucose uptake assay to them is a nonsensical procedure that will not yield any meaningful data about mitochondrial function.

**D) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
*   **Principle:** The long chemical name refers to the common and well-established fluorescent dye **JC-1**. JC-1 is a cationic dye used to measure the mitochondrial membrane potential (ΔΨm).
*   **Relevance:** The mitochondrial membrane potential is a critical indicator of mitochondrial health. It is generated by the electron transport chain and is essential for ATP synthesis. In healthy mitochondria with a high potential, JC-1 forms red fluorescent aggregates. If a drug causes dysfunction and the potential collapses, the dye disperses as green fluorescent monomers. Flow cytometry allows for a rapid and quantitative measurement of this red-to-green shift across thousands of cells.
*   **Careful Point:** This is a standard, powerful, and direct method for assessing a drug's impact on mitochondrial integrity and function. It is a very helpful experiment.

### Conclusion

Experiments A, B, and D all describe valid methods to measure key parameters directly or indirectly related to mitochondrial function: ATP levels/cell viability, specific intramitochondrial activity, and membrane potential, respectively. Experiment C is the only one that is fundamentally flawed because it attempts to measure a process (glucose uptake) on an organelle (the mitochondrion) that does not perform that function.

<<<C>>>
"""

# Run the check
result = check_mitochondria_experiment_answer(llm_answer)
print(result)