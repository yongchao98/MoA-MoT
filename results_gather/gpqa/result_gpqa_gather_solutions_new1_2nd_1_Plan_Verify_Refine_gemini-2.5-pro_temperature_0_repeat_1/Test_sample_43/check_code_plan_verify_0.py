import re

def check_answer(question, candidate_answers):
    """
    Checks the correctness of the answer to a biological question about mitochondrial experiments.

    The function encodes key biological facts about mitochondrial metabolism and evaluates
    each experimental option against these facts.
    """

    # --- Knowledge Base ---
    # Define the core biological principles relevant to the question.
    class BiologicalFacts:
        MITOCHONDRIAL_FUNCTIONS = {
            "ATP Production",
            "Membrane Potential Maintenance",
            "Metabolic Activity (ROS/Heat)",
        }
        # Mitochondria import pyruvate, the product of glycolysis. They do NOT import glucose.
        MITOCHONDRIAL_IMPORTS = {"pyruvate", "fatty acids"}
        PLASMA_MEMBRANE_PROCESSES = {"glucose uptake"}
        CYTOPLASMIC_PROCESSES = {"glycolysis (glucose -> pyruvate)"}

    # --- Experimental Options ---
    # Represent each experimental option with its key characteristics.
    # The 'reagent' field helps identify the specific technique.
    experiments = {
        'A': {
            "name": "Luciferase Assay",
            "measures": "ATP",
            "reagent": "luciferase/luciferin",
            "is_helpful": True,
            "reasoning": "Measures ATP, the primary output of mitochondrial energy production. It is a valid, though sometimes indirect, indicator of mitochondrial function."
        },
        'B': {
            "name": "Glucose Uptake on Isolated Mitochondria",
            "measures": "Glucose Uptake",
            "reagent": "Glucose Uptake Colorimetric Assay Kit",
            "is_helpful": False, # This is the flawed experiment
            "reasoning": f"This experiment is fundamentally flawed. 'Glucose uptake' is a process of the plasma membrane, not mitochondria. Mitochondria import {', '.join(BiologicalFacts.MITOCHONDRIAL_IMPORTS)}, not glucose."
        },
        'C': {
            "name": "Mitochondrial Membrane Potential Assay",
            "measures": "Mitochondrial Membrane Potential",
            "reagent": "5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide (JC-1)",
            "is_helpful": True,
            "reasoning": "Measures mitochondrial membrane potential, a direct and critical indicator of mitochondrial health and function."
        },
        'D': {
            "name": "Mitochondrial Activity/Stress Imaging",
            "measures": "Mitochondrial-specific parameters (e.g., ROS, temperature)",
            "reagent": "Mito-RTP stain",
            "is_helpful": True,
            "reasoning": "Uses a mitochondria-targeted probe to directly measure parameters related to metabolic activity or stress."
        }
    }

    # --- Evaluation ---
    # Find the experiment that is NOT helpful based on our knowledge base.
    unhelpful_experiment_key = None
    unhelpful_experiment_reason = ""
    for key, exp in experiments.items():
        if not exp["is_helpful"]:
            unhelpful_experiment_key = key
            unhelpful_experiment_reason = exp["reasoning"]
            break

    # --- Check the provided answer ---
    # Extract the final letter from the last candidate answer.
    # The provided answer block has the format <<<B>>>
    try:
        # Find the last occurrence of <<<...>>>
        last_answer_match = re.findall(r'<<<([A-Z])>>>', candidate_answers)
        if not last_answer_match:
            return "Could not find a final answer in the format <<<A>>> in the provided text."
        
        provided_answer_key = last_answer_match[-1]

        if provided_answer_key not in experiments:
            return f"The provided answer '{provided_answer_key}' is not a valid option (A, B, C, or D)."

        # Compare the provided answer with the one identified by our logic.
        if provided_answer_key == unhelpful_experiment_key:
            return "Correct"
        else:
            return (f"Incorrect. The provided answer is {provided_answer_key}, but the correct answer is {unhelpful_experiment_key}. "
                    f"Reason: {unhelpful_experiment_reason}")

    except Exception as e:
        return f"An error occurred during parsing or checking: {e}"

# The question and candidate answers are provided as multiline strings.
question = """
A eukaryotic cell evolved a mechanism to turn macromolecular building blocks into energy. The process occurs in mitochondria, which are cellular energy factories. In the series of redox reactions, the energy from food is stored between the phosphate groups and used as a universal cellular currency. The energy-laden molecules are shuttled out of the mitochondrion to serve in all cellular processes. You discovered a new anti-diabetes drug and want to investigate whether it has an effect on the mitochondria. You set up a range of experiments with your HEK293 cell line. Which of the experiments listed below will not help you discover the mitochondrial role of your drug:

A) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant
B) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit
C) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide
D) Confocal fluorescence microscopy after Mito-RTP staining of the cells
"""

candidate_answers = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

**1. Synthesize the Consensus from Candidate Answers**

All 14 candidate answers, despite sometimes using different lettering for the options, arrive at the same conclusion through identical biological reasoning. The unanimous consensus is that the experiment designed to measure glucose uptake in isolated mitochondria is fundamentally flawed and will not provide useful information. The core reasoning is that mitochondria do not perform glucose uptake; this process occurs at the cell's plasma membrane, and mitochondria import the downstream product, pyruvate.

**2. Analyze Each Option from the Original Question**

Let's evaluate each choice based on the robust reasoning provided by the candidates.

*   **A) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
    *   **What it measures:** This is an assay for ATP. The luciferase enzyme requires ATP to produce light.
    *   **Careful Point:** Mitochondria are the primary producers of cellular ATP. Measuring ATP levels is a fundamental way to assess their energy-producing function. While measuring ATP in the supernatant is an indirect method, it can serve as a valid indicator of cell death caused by severe mitochondrial toxicity.
    *   **Conclusion:** This experiment is helpful for investigating the drug's effect.

*   **B) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
    *   **What it measures:** This experiment first isolates mitochondria and then attempts to measure their ability to take up glucose.
    *   **Careful Point (The Critical Flaw):** This experiment is based on a flawed biological premise. Glucose is transported across the cell's outer membrane and broken down into pyruvate in the cytoplasm (glycolysis). It is **pyruvate**, not glucose, that is then transported into the mitochondria. Mitochondria lack the necessary transporters to take up glucose directly.
    *   **Conclusion:** This experiment is designed to measure a process that does not occur in the target organelle. It is biologically nonsensical and **will not be helpful**.

*   **C) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
    *   **What it measures:** The chemical name is for the common dye **JC-1**, which is used to measure the mitochondrial membrane potential (ΔΨm).
    *   **Careful Point:** The mitochondrial membrane potential is a direct and critical indicator of the health of the electron transport chain and the mitochondrion's ability to produce ATP. A drop in this potential is a classic sign of mitochondrial dysfunction.
    *   **Conclusion:** This is a standard, powerful, and highly relevant experiment that is helpful.

*   **D) Confocal fluorescence microscopy after Mito-RTP staining of the cells**
    *   **What it measures:** The "Mito-" prefix indicates the stain specifically targets mitochondria. Probes like this are used to measure key parameters like reactive oxygen species (ROS), temperature (a proxy for metabolic rate), or morphology.
    *   **Careful Point:** Measuring these parameters provides a direct look at the functional state and stress levels of the mitochondria within an intact cell.
    *   **Conclusion:** This is a direct and valid method for investigation and is helpful.

**Final Decision**

Experiments A, C, and D all measure valid parameters related to mitochondrial function (ATP production, membrane potential, and metabolic activity/stress). Experiment B is the only one based on a fundamentally incorrect biological premise and will therefore not yield any useful information.

<<<B>>>
"""

# Run the check
result = check_answer(question, candidate_answers)
print(result)