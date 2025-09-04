import re

def check_answer_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given biological question.

    The function encodes the biological principles behind each experimental option
    to determine which one is invalid for studying mitochondrial function.

    Args:
        question (str): The question text (not used in this specific implementation but good practice).
        llm_answer_text (str): The full text of the LLM's response, including the final answer in <<<>>> format.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """

    # --- Knowledge Base ---
    # This dictionary encodes the scientific validity of each experimental option.
    # The question asks which experiment will NOT help.
    experiments = {
        'A': {
            'description': "Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit",
            'is_helpful': False,
            'reasoning': "This experiment is fundamentally flawed. Mitochondria do not have transporters to take up glucose directly. They import pyruvate, the end-product of glycolysis which occurs in the cytoplasm. Therefore, a glucose uptake assay on isolated mitochondria is biologically nonsensical."
        },
        'B': {
            'description': "Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide",
            'is_helpful': True,
            'reasoning': "This describes using the JC-1 dye to measure mitochondrial membrane potential (ΔΨm). This is a standard, direct, and quantitative method to assess mitochondrial health and function."
        },
        'C': {
            'description': "Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant",
            'is_helpful': True,
            'reasoning': "This is a luciferase assay to measure ATP. Since mitochondria are the primary producers of ATP, this is a relevant measure of their function. Measuring extracellular ATP can also indicate cell death caused by mitochondrial toxins."
        },
        'D': {
            'description': "Confocal fluorescence microscopy after Mito-RTP staining of the cells",
            'is_helpful': True,
            'reasoning': "This uses a mitochondria-specific ('Mito-') probe to measure parameters like ROS, temperature, or morphology. This is a direct way to visualize and quantify a drug's effect on mitochondrial activity or stress."
        }
    }

    # --- Logic to find the correct answer ---
    # The correct answer is the key of the experiment that is NOT helpful.
    correct_key = None
    for key, data in experiments.items():
        if not data['is_helpful']:
            correct_key = key
            break

    if correct_key is None:
        return "Error in checker logic: No single unhelpful experiment was identified."

    # --- Extract the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The answer format is invalid. It should be '<<<X>>>' where X is A, B, C, or D."

    llm_key = match.group(1)

    # --- Compare and generate result ---
    if llm_key == correct_key:
        return "Correct"
    else:
        reason_for_llm_choice_being_wrong = experiments[llm_key]['reasoning']
        reason_for_correct_choice = experiments[correct_key]['reasoning']
        
        return (f"Incorrect. The provided answer is {llm_key}, but the correct answer is {correct_key}. "
                f"The question asks which experiment will NOT help. "
                f"Option {llm_key} is actually a helpful experiment because: {reason_for_llm_choice_being_wrong}. "
                f"The correct answer is {correct_key} because: {reason_for_correct_choice}")


# Provided question and answer for the check
question = """
A eukaryotic cell evolved a mechanism to turn macromolecular building blocks into energy. The process occurs in mitochondria, which are cellular energy factories. In the series of redox reactions, the energy from food is stored between the phosphate groups and used as a universal cellular currency. The energy-laden molecules are shuttled out of the mitochondrion to serve in all cellular processes. You discovered a new anti-diabetes drug and want to investigate whether it has an effect on the mitochondria. You set up a range of experiments with your HEK293 cell line. Which of the experiments listed below will not help you discover the mitochondrial role of your drug:

A) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit
B) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide
C) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant
D) Confocal fluorescence microscopy after Mito-RTP staining of the cells
"""

llm_answer = """
Here is a step-by-step analysis of the experimental options to determine the final answer. The goal is to identify which experiment will **not** be helpful for investigating a drug's effect on mitochondria.

**Step 1: Analyze Option A - Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
*   **Procedure:** This experiment first isolates mitochondria from the rest of the cell's components. Then, it attempts to measure the rate at which this isolated fraction takes up glucose.
*   **Biological Principle:** The process of converting glucose to energy begins with glycolysis, where glucose is broken down into pyruvate. This entire process occurs in the cell's cytoplasm. The initial step, glucose uptake from the environment, is performed by glucose transporters (e.g., GLUTs) located on the cell's outer plasma membrane.
*   **Careful Point:** Mitochondria do not have the necessary transporters to take up glucose directly. They import the end-product of glycolysis, pyruvate, from the cytoplasm. Therefore, applying a glucose uptake assay to isolated mitochondria is based on a fundamentally incorrect biological premise.
*   **Conclusion:** This experiment is designed to measure a process in an organelle that does not perform it. The results would be meaningless, so this experiment **will not be helpful**.

**Step 2: Analyze Option B - Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
*   **Procedure:** This uses a specific chemical dye followed by analysis with a flow cytometer.
*   **Biological Principle:** The long chemical name corresponds to the well-known fluorescent dye **JC-1**. JC-1 is a cationic dye used specifically to measure the mitochondrial membrane potential (ΔΨm), a key indicator of mitochondrial health. In healthy mitochondria, the dye forms aggregates that fluoresce red; in unhealthy mitochondria with a collapsed potential, it remains as monomers that fluoresce green.
*   **Careful Point:** Flow cytometry allows for the rapid and quantitative measurement of the red-to-green fluorescence ratio in thousands of individual cells. This provides a robust and direct assessment of the drug's impact on mitochondrial function.
*   **Conclusion:** This is a standard, direct, and powerful method for investigating a drug's effect on mitochondria. This experiment **will be helpful**.

**Step 3: Analyze Option C - Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
*   **Procedure:** This describes a luciferase-based assay to measure the concentration of ATP. The enzyme luciferase uses ATP to produce light in the presence of its substrate, luciferin.
*   **Biological Principle:** Mitochondria are the primary producers of ATP in the cell. A drug that impairs mitochondrial function would almost certainly lead to a change in the cell's overall energy status and ATP levels.
*   **Careful Point:** This experiment measures ATP in the supernatant (the extracellular space). While a direct measurement would lyse the cells and measure intracellular ATP, a significant release of ATP into the supernatant is often a sign of cellular stress or death (necrosis), which can be triggered by mitochondrial toxins. Thus, it can serve as an indirect, but still potentially useful, indicator of mitochondrial toxicity.
*   **Conclusion:** This experiment measures a parameter (ATP) that is directly related to the primary function of mitochondria and **will be helpful**.

**Step 4: Analyze Option D - Confocal fluorescence microscopy after Mito-RTP staining of the cells**
*   **Procedure:** This involves staining cells with a mitochondria-specific probe ("Mito-RTP") and visualizing them with a high-resolution confocal microscope.
*   **Biological Principle:** Probes with a "Mito-" prefix are designed to specifically accumulate in mitochondria. These probes can be used to assess various parameters, such as mitochondrial mass, morphology, temperature, or the production of reactive oxygen species (ROS). All of these are directly related to mitochondrial health and function.
*   **Careful Point:** Confocal microscopy allows for detailed visualization of the mitochondrial network and the quantification of the probe's signal within the organelle, providing a direct assessment of how the drug affects mitochondrial structure and/or activity.
*   **Conclusion:** This is a valid and direct method for investigating the mitochondrial role of a drug. This experiment **will be helpful**.

**Final Decision:**
Comparing the four options, experiments B, C, and D all measure valid parameters of mitochondrial function or health (membrane potential, ATP levels/cell death, and morphology/activity, respectively). Experiment A is the only one based on a flawed biological premise, as it attempts to measure glucose uptake in an organelle that does not perform this function.

<<<A>>>
"""

# Run the check
result = check_answer_correctness(question, llm_answer)
print(result)