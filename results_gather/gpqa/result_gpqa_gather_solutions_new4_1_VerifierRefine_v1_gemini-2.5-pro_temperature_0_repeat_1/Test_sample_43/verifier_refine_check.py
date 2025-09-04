import re

def check_answer_correctness(question_text, llm_answer_text):
    """
    Checks the correctness of an LLM's answer to a multiple-choice question about mitochondrial experiments.

    Args:
        question_text (str): The text of the question with the options.
        llm_answer_text (str): The text of the LLM's answer, including its reasoning and final choice.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """

    # Step 1: Define the ground truth based on established cell biology principles.
    # The question asks to identify the experiment that will NOT help.
    
    # Fact A: JC-1 measures mitochondrial membrane potential, a key indicator of mitochondrial health. This is a HELPFUL experiment.
    # Fact B: Mitochondria do NOT directly take up glucose. Glucose is converted to pyruvate in the cytoplasm, and pyruvate enters the mitochondria.
    #         Therefore, a glucose uptake assay on isolated mitochondria is biologically incorrect. This is NOT a helpful experiment.
    # Fact C: "Mito-" prefixed probes (like Mito-RTP) are designed to assess mitochondrial parameters (e.g., temperature, ROS, morphology). This is a HELPFUL experiment.
    # Fact D: Luciferase assays measure ATP. Since mitochondria are the primary producers of ATP, this is a fundamental way to assess their function. This is a HELPFUL experiment.

    correct_option = 'B'
    reasoning = {
        'A': "This experiment is helpful. The dye JC-1 is a standard probe for measuring mitochondrial membrane potential, which is a key indicator of mitochondrial health and function.",
        'B': "This experiment is NOT helpful. The experimental design is fundamentally flawed. Mitochondria do not directly take up glucose; this process occurs at the cell's plasma membrane. The product of cytoplasmic glucose breakdown (pyruvate) is what enters the mitochondria. Therefore, a glucose uptake assay on isolated mitochondria will not yield meaningful data.",
        'C': "This experiment is helpful. 'Mito-RTP' is a mitochondria-targeted probe used to measure a parameter related to metabolic activity (like temperature or reactive oxygen species). This is a valid method to assess mitochondrial function.",
        'D': "This experiment is helpful. The luciferase assay measures ATP levels. Since mitochondria are the primary producers of cellular ATP, measuring ATP is a direct way to assess the impact of a drug on mitochondrial energy production."
    }

    # Step 2: Extract the final answer from the LLM's response.
    # The expected format is <<<X>>>.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."

    llm_choice = match.group(1)

    # Step 3: Compare the LLM's answer with the ground truth.
    if llm_choice == correct_option:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer is <<<{llm_choice}>>>.\n"
            f"The correct answer is <<<{correct_option}>>>.\n\n"
            f"Reasoning:\n"
            f"- The chosen answer <<<{llm_choice}>>> is incorrect because: {reasoning[llm_choice]}\n"
            f"- The correct answer <<<{correct_option}>>> is correct because: {reasoning[correct_option]}"
        )
        return error_message

# The user's provided question and answer text would be passed to this function.
# For this example, I will use the provided inputs.

question = """
A eukaryotic cell evolved a mechanism to turn macromolecular building blocks into energy. The process occurs in mitochondria, which are cellular energy factories. In the series of redox reactions, the energy from food is stored between the phosphate groups and used as a universal cellular currency. The energy-laden molecules are shuttled out of the mitochondrion to serve in all cellular processes. You discovered a new anti-diabetes drug and want to investigate whether it has an effect on the mitochondria. You set up a range of experiments with your HEK293 cell line. Which of the experiments listed below will not help you discover the mitochondrial role of your drug:

A) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide
B) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit
C) Confocal fluorescence microscopy after Mito-RTP staining of the cells
D) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant
"""

llm_answer = """
Here is a step-by-step analysis of the question and the provided options to determine the final answer.

**1. Deconstruct the Question**
The core task is to identify which of the four proposed experiments is **not** suitable for investigating a drug's effect on mitochondria. The context is cellular energy production, where mitochondria convert building blocks into ATP.

**2. Analyze Each Experimental Option**

*   **Option A: Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
    *   **Reagent:** The long chemical name refers to the common fluorescent dye **JC-1**.
    *   **Function:** JC-1 is a potentiometric dye used to measure the mitochondrial membrane potential (ΔΨm). A high membrane potential is a key indicator of healthy, functioning mitochondria. In healthy mitochondria, JC-1 forms aggregates that fluoresce red; in unhealthy mitochondria with a low potential, it remains as monomers that fluoresce green.
    *   **Technique:** Flow cytometry quantifies this red/green fluorescence ratio in thousands of individual cells.
    *   **Conclusion:** This is a direct, standard, and powerful method to assess mitochondrial health. This experiment **is helpful**.

*   **Option B: Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
    *   **Step 1: Differential centrifugation.** This is a valid technique to isolate mitochondria from other cellular components.
    *   **Step 2: Glucose Uptake Assay.** This assay measures the transport of glucose into a biological system.
    *   **Biological Flaw:** The critical point is that glucose uptake occurs at the cell's plasma membrane, not the mitochondrial membrane. Glucose is broken down into pyruvate in the cytoplasm during glycolysis. It is the **pyruvate**, not glucose, that is then transported into the mitochondria. Isolated mitochondria lack the necessary transporters (like GLUTs) to take up glucose.
    *   **Conclusion:** This experimental design is fundamentally flawed. Applying a glucose uptake assay to isolated mitochondria is biochemically incorrect and will not yield any meaningful information about mitochondrial function. This experiment **is not helpful**.

*   **Option C: Confocal fluorescence microscopy after Mito-RTP staining of the cells**
    *   **Reagent:** "Mito-RTP" is a mitochondria-targeted ("Mito-") probe. Research indicates it can be a probe for reactive oxygen species (ROS) or temperature, both of which are directly related to mitochondrial metabolic activity.
    *   **Function:** The electron transport chain in mitochondria is a major source of ROS, and its activity generates heat. Measuring changes in either parameter provides insight into mitochondrial function and stress.
    *   **Technique:** Confocal microscopy provides high-resolution images, allowing for the visualization and quantification of the probe's signal specifically within the mitochondria.
    *   **Conclusion:** This is a direct and valid method to investigate the drug's effect on mitochondrial activity. This experiment **is helpful**.

*   **Option D: Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
    *   **Reagent/Principle:** This is a luciferase-based ATP assay. The enzyme luciferase uses ATP to produce light. The amount of light is proportional to the amount of ATP.
    *   **Function:** Mitochondria are the primary producers of cellular ATP. Measuring ATP levels is a fundamental way to assess the primary output of mitochondrial function. A drug that impairs mitochondria would almost certainly lead to a change in cellular ATP levels.
    *   **Conclusion:** This is a highly relevant method for assessing the drug's impact on cellular energy metabolism, which is governed by the mitochondria. This experiment **is helpful**.

**3. Synthesize and Decide**
Experiments A, C, and D all measure key parameters directly or indirectly related to mitochondrial health and function: membrane potential, metabolic activity/stress, and ATP production, respectively.

Experiment B is the only one based on a flawed biological premise. It attempts to measure a process (glucose uptake) on an organelle (the mitochondrion) that does not perform that function. Therefore, it is the only experiment that will not help discover the mitochondrial role of the drug.

<<<B>>>
"""

# Run the check
result = check_answer_correctness(question, llm_answer)
print(result)