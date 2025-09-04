import re

def check_answer(llm_output):
    """
    Checks the correctness of the LLM's answer to the biology question.

    This function encodes the biological principles of cellular respiration and mitochondrial function
    to logically determine which experiment is not helpful. It then compares this
    conclusion to the LLM's provided answer.
    """

    # --- Knowledge Base ---
    # Define the known functions of mitochondria and other cellular components.
    knowledge_base = {
        "mitochondrial_functions": [
            "ATP_production",  # The primary energy currency production.
            "membrane_potential",  # Essential for ATP synthesis, measured by JC-1.
            "metabolic_activity",  # Generates heat and ROS, measured by probes like Mito-RTP.
            "pyruvate_import"  # Imports the product of glycolysis.
        ],
        "plasma_membrane_functions": [
            "glucose_uptake"  # Mediated by GLUT transporters on the cell surface.
        ],
        "cytoplasm_functions": [
            "glycolysis"  # Breakdown of glucose into pyruvate.
        ]
    }

    # --- Experimental Analysis ---
    # Map each experimental option to what it measures and its target.
    experiments = {
        "A": {
            "description": "Confocal fluorescence microscopy after Mito-RTP staining",
            "measures": "metabolic_activity",  # Measures heat/ROS in mitochondria.
            "target": "mitochondria",
            "is_helpful": True  # This is a direct measure of a mitochondrial function.
        },
        "B": {
            "description": "Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit",
            "measures": "glucose_uptake",
            "target": "isolated_mitochondria",
            # This is the key check: Is the function measured performed by the target?
            "is_helpful": "glucose_uptake" not in knowledge_base["mitochondrial_functions"] and "glucose_uptake" in knowledge_base["plasma_membrane_functions"]
        },
        "C": {
            "description": "Transformation of cells with recombinant luciferase and luminometer reading",
            "measures": "ATP_production",  # Measures the primary output of mitochondria.
            "target": "cell",
            "is_helpful": True # This is a direct measure of a key mitochondrial output.
        },
        "D": {
            "description": "Flow cytometry after labeling with JC-1 dye",
            "measures": "membrane_potential", # JC-1 measures mitochondrial membrane potential.
            "target": "mitochondria",
            "is_helpful": True # This is a direct measure of mitochondrial health.
        }
    }

    # --- Logic to find the unhelpful experiment ---
    unhelpful_experiment_key = None
    reason = ""
    for key, details in experiments.items():
        # The logic for B is special: it's helpful if the function is mitochondrial.
        # Since glucose_uptake is NOT a mitochondrial function, is_helpful becomes False.
        if not details["is_helpful"]:
            unhelpful_experiment_key = key
            reason = (f"Experiment {key} ('{details['description']}') is not helpful because it "
                      f"attempts to measure '{details['measures']}' on '{details['target']}'. "
                      f"However, '{details['measures']}' is a function of the plasma membrane, "
                      f"not the mitochondria. Mitochondria import pyruvate, not glucose.")
            break

    # --- Final Check ---
    # Extract the answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_output)
    if not match:
        return "Error: Could not find the answer in the required format '<<<...>>>' in the provided text."

    provided_answer = match.group(1)

    if provided_answer == unhelpful_experiment_key:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is <<<{provided_answer}>>>, but this experiment is a valid method to study mitochondria.\n"
                f"The correct answer is <<<{unhelpful_experiment_key}>>>.\n"
                f"Reason: {reason}")


# The final response from the LLM being checked
llm_response = """
Here is a step-by-step analysis of the question to determine the final answer.

**1. Understand the Goal of the Question**
The question asks to identify which of the four listed experimental methods is **not** suitable for investigating a drug's effect on mitochondria. This means we are looking for an experiment that is either irrelevant to mitochondrial function or is based on a flawed biological premise.

**2. Analyze Each Experimental Option**

*   **A) Confocal fluorescence microscopy after Mito-RTP staining of the cells**
    *   **Principle:** This experiment uses a fluorescent probe with the prefix "Mito-", indicating it is designed to specifically accumulate within mitochondria. Probes like this (e.g., Mito-RTP, MitoTracker, MitoSOX) are used to measure specific parameters of mitochondrial function or status, such as temperature (a proxy for metabolic rate), reactive oxygen species (ROS), or membrane potential.
    *   **Careful Point:** Measuring changes in metabolic activity (heat) or oxidative stress (ROS) directly within the mitochondria is a direct and valid way to assess a drug's impact.
    *   **Conclusion:** This experiment is helpful.

*   **B) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit**
    *   **Principle:** This experiment involves two steps. First, mitochondria are physically isolated from the rest of the cell. Second, a glucose uptake assay is performed on this isolated mitochondrial fraction.
    *   **Careful Point (The Flaw):** This experimental design is based on a fundamental biological error. The process of glucose uptake from the environment is carried out by glucose transporters (like GLUTs) located on the **cell's plasma membrane**, not on the mitochondrial membranes. The initial breakdown of glucose (glycolysis) occurs in the **cytoplasm**, producing pyruvate. It is **pyruvate**, not glucose, that is then transported into the mitochondria to fuel the Krebs cycle.
    *   **Conclusion:** Applying a glucose uptake assay to isolated mitochondria is a nonsensical procedure. It will not yield any meaningful information about mitochondrial function. This experiment is **not helpful**.

*   **C) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant**
    *   **Principle:** This describes a luciferase-based assay, which measures the concentration of ATP. The enzyme luciferase uses ATP to produce light from its substrate, luciferin.
    *   **Careful Point:** Mitochondria are the primary producers of ATP in the cell via oxidative phosphorylation. Therefore, measuring the level of ATP is a direct assessment of the primary output of mitochondrial function. A drug that impairs mitochondria would almost certainly lead to a change in cellular ATP levels. Even if measuring extracellular ATP (in the supernatant), it can be an indicator of cell stress or death caused by mitochondrial damage.
    *   **Conclusion:** This experiment is helpful.

*   **D) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide**
    *   **Principle:** The long chemical name refers to the common fluorescent dye **JC-1**. JC-1 is a probe used to measure the mitochondrial membrane potential (ΔΨm), a critical indicator of mitochondrial health. In healthy, energized mitochondria, JC-1 forms red fluorescent aggregates. If the potential collapses due to dysfunction, the dye disperses as green fluorescent monomers.
    *   **Careful Point:** The mitochondrial membrane potential is directly generated by the redox reactions of the electron transport chain and is essential for ATP synthesis. This is a standard and powerful method.
    *   **Conclusion:** This experiment is helpful.

**3. Final Decision**
Experiments A, C, and D all measure valid parameters directly or indirectly related to mitochondrial function (metabolic activity, ATP production, and membrane potential, respectively). Experiment B is based on a flawed biological premise, as isolated mitochondria do not perform glucose uptake. Therefore, it is the only experiment that will not provide useful information.

<<<B>>>
"""

# Run the check
result = check_answer(llm_response)
print(result)