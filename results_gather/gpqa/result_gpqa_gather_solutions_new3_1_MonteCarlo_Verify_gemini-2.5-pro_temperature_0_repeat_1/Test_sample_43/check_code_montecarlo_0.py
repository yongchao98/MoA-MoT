import re

def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer for the given cell biology question.

    The function defines the ground truth for each experimental option based on
    established biological principles and compares it with the provided answer.
    """

    # The final answer provided by the LLM.
    llm_answer_text = """
    Here is a step-by-step analysis of the question and the experimental options provided.

    **Objective**: Identify the single experiment from the list that is **not** suitable for investigating a drug's effect on mitochondrial function.

    **Step 1: Analyze the biological context.**
    The question describes the role of mitochondria in converting "macromolecular building blocks into energy" via redox reactions, storing this energy in ATP. This refers to cellular respiration, specifically the Krebs cycle and oxidative phosphorylation, which occur inside the mitochondria. The key inputs for this process are not glucose itself, but rather pyruvate (from glycolysis in the cytoplasm), fatty acids, and amino acids. The key outputs and indicators of function are ATP production, maintenance of a high mitochondrial membrane potential (ΔΨm), and controlled levels of byproducts like reactive oxygen species (ROS) and heat.

    **Step 2: Evaluate each experimental option.**

    *   **A) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide:**
        *   **What it is:** The long chemical name refers to the common fluorescent dye **JC-1**.
        *   **What it measures:** JC-1 is a cationic dye used to measure mitochondrial membrane potential (ΔΨm). In healthy, high-potential mitochondria, it forms red fluorescent aggregates. In unhealthy, low-potential mitochondria, it remains as green fluorescent monomers.
        *   **Careful Point:** The mitochondrial membrane potential is a direct and critical indicator of mitochondrial health and its ability to produce ATP.
        *   **Conclusion:** This is a standard, powerful, and highly relevant experiment. It **would be helpful**.

    *   **B) Confocal fluorescence microscopy after Mito-RTP staining of the cells:**
        *   **What it is:** An imaging experiment using a mitochondria-targeted ("Mito-") fluorescent probe. Probes like this are designed to report on specific conditions within the mitochondria.
        *   **What it measures:** Depending on the specific probe (e.g., MitoSOX for ROS, Mito-Thermo-Yellow for temperature), this technique can measure key indicators of mitochondrial activity or stress.
        *   **Careful Point:** Visualizing changes in mitochondrial morphology, stress levels (like ROS), or metabolic activity (like temperature) provides direct evidence of a drug's effect at the subcellular level.
        *   **Conclusion:** This is a valid and direct method for investigation. It **would be helpful**.

    *   **C) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit:**
        *   **What it is:** A two-step experiment. First, mitochondria are isolated. Second, a glucose uptake assay is performed on this isolated fraction.
        *   **What it measures:** The assay measures the transport of glucose into a biological compartment.
        *   **Careful Point:** This experiment has a fundamental design flaw. Glucose uptake is a process that occurs at the cell's plasma membrane, not the mitochondrial membrane. The cell takes up glucose into the cytoplasm, where it is converted to pyruvate via glycolysis. It is the pyruvate that is then transported into the mitochondria. Isolated mitochondria do not have the necessary transporters to take up glucose.
        *   **Conclusion:** This experiment is based on an incorrect biological premise and will not yield any meaningful information about mitochondrial function. It **will not be helpful**.

    *   **D) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant:**
        *   **What it is:** A luciferase-based assay to measure ATP. The luciferase enzyme uses ATP to produce light from luciferin.
        *   **What it measures:** This assay quantifies ATP levels.
        *   **Careful Point:** Since mitochondria are the primary producers of cellular ATP, measuring ATP concentration is a fundamental readout of the cell's overall energy status and, by extension, mitochondrial function. A drug that impairs mitochondria would almost certainly affect cellular ATP levels.
        *   **Conclusion:** This is a highly relevant experiment for assessing the primary output of mitochondria. It **would be helpful**.

    **Final Decision:**
    Experiments A, B, and D all measure valid parameters directly or indirectly related to mitochondrial function (membrane potential, metabolic activity/stress, and ATP production). Experiment C attempts to measure a process (glucose uptake) in an organelle (the mitochondrion) that does not perform it, making the experimental design fundamentally flawed and unhelpful for the stated goal.

    <<<C>>>
    """

    # Ground truth based on biological principles. The question asks which experiment will NOT help.
    experiments_validity = {
        'A': {
            'is_helpful': True,
            'reason': "It measures mitochondrial membrane potential, a key indicator of mitochondrial health."
        },
        'B': {
            'is_helpful': True,
            'reason': "It uses a mitochondria-targeted probe to assess mitochondrial activity or stress."
        },
        'C': {
            'is_helpful': False,
            'reason': "It is based on a flawed premise. Mitochondria do not take up glucose; this process occurs at the cell's plasma membrane and in the cytoplasm."
        },
        'D': {
            'is_helpful': True,
            'reason': "It measures ATP levels, the primary energy product of mitochondria."
        }
    }

    # The correct answer is the option that is NOT helpful.
    correct_option = [opt for opt, details in experiments_validity.items() if not details['is_helpful']][0]

    # Parse the LLM's final answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<X>>>'."

    llm_choice = match.group(1)

    # Check if the LLM's choice matches the correct option
    if llm_choice == correct_option:
        return "Correct"
    else:
        reasoning = (f"Incorrect: The answer '{llm_choice}' is wrong. "
                     f"The question asks for the experiment that will NOT help investigate mitochondrial function.\n"
                     f"Experiment '{llm_choice}' is actually a valid method because {experiments_validity[llm_choice]['reason']}\n"
                     f"The correct answer is '{correct_option}' because {experiments_validity[correct_option]['reason']}")
        return reasoning

# Execute the check and print the result
print(check_answer_correctness())