import re

def check_answer_correctness(llm_answer_text):
    """
    Checks the correctness of the LLM's answer by applying biological principles.

    The function evaluates four experimental options to determine which one is NOT suitable
    for studying mitochondrial function. It then compares this correct choice with the
    provided LLM answer.

    Args:
        llm_answer_text (str): The full response from the LLM, including the final answer in <<<>>> format.

    Returns:
        str: "Correct" if the answer is right, or a string explaining the error.
    """

    # Extract the single-letter answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<A>>>."
    llm_answer = match.group(1)

    # Database of biological facts for each experimental option
    experiments = [
        {
            "id": "A",
            "description": "Flow cytometry after labeling with JC-1 dye",
            "measures": "Mitochondrial Membrane Potential (ΔΨm)",
            "is_valid_for_mitochondria": True,
            "reasoning": "JC-1 is a standard probe for mitochondrial membrane potential, a key indicator of mitochondrial health and function."
        },
        {
            "id": "B",
            "description": "Confocal fluorescence microscopy after Mito-RTP staining",
            "measures": "Mitochondrial-specific parameters (e.g., temperature, ROS, pyrophosphate)",
            "is_valid_for_mitochondria": True,
            "reasoning": "Mito-probes are designed to measure parameters like metabolic activity or stress directly within the mitochondria."
        },
        {
            "id": "C",
            "description": "Luciferase assay for ATP",
            "measures": "ATP levels",
            "is_valid_for_mitochondria": True,
            "reasoning": "Mitochondria are the primary producers of cellular ATP, so measuring ATP is a fundamental way to assess their energy output."
        },
        {
            "id": "D",
            "description": "Glucose Uptake Colorimetric Assay on isolated mitochondria",
            "measures": "Glucose Uptake",
            "is_valid_for_mitochondria": False,
            "reasoning": "This experiment is fundamentally flawed. Glucose uptake occurs at the cell's plasma membrane, not in the mitochondria. Mitochondria import pyruvate (the product of glycolysis), not glucose. Therefore, this assay will not provide useful information about mitochondrial function."
        }
    ]

    # Identify the correct answer based on biological principles
    # The correct answer to the question is the experiment that is NOT valid.
    correct_option = None
    invalid_reason = ""
    for exp in experiments:
        if not exp["is_valid_for_mitochondria"]:
            correct_option = exp["id"]
            invalid_reason = exp["reasoning"]
            break

    # Check if the LLM's answer matches the correct option
    if llm_answer == correct_option:
        return "Correct"
    else:
        # Find the reason why the LLM's chosen answer is actually a valid experiment
        llm_choice_info = next((exp for exp in experiments if exp["id"] == llm_answer), None)
        if llm_choice_info:
            llm_choice_reason = llm_choice_info["reasoning"]
            return (f"Incorrect. The provided answer '{llm_answer}' is wrong. "
                    f"The experiment in option '{llm_answer}' is a valid method to study mitochondria because {llm_choice_reason}. "
                    f"The correct answer is '{correct_option}' because the experiment described is biologically flawed: {invalid_reason}")
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

# The final answer provided by the LLM to be checked
llm_final_answer = """
Here is a step-by-step analysis to decide the final answer.

**1. Deconstruct the Question**
The core task is to identify which of the four proposed experiments is **not suitable** for investigating a drug's effect on mitochondria. The background information establishes that mitochondria are the cell's "energy factories," using redox reactions to produce energy stored in ATP. This means any valid experiment should measure a parameter related to this process.

**2. Analyze Each Experimental Option**

*   **A) Flow cytometry after labeling with 2.5 μM 5,5',6,6'-Tetrachloro-1,1',3, 3'-tetraethylbenzimidazolylcarbocyanine iodide:**
    *   **Careful Point 1: Identify the reagent.** The long chemical name is for the common fluorescent dye **JC-1**.
    *   **Careful Point 2: Understand the mechanism.** JC-1 is a cationic dye used specifically to measure mitochondrial membrane potential (ΔΨm). This potential is a direct result of the electron transport chain's activity and is essential for ATP synthesis. In healthy, energized mitochondria, JC-1 forms red fluorescent aggregates. In unhealthy, depolarized mitochondria, it remains as green fluorescent monomers.
    *   **Conclusion:** Measuring the red/green fluorescence ratio with flow cytometry is a direct, standard, and powerful method to assess mitochondrial health. This experiment **is helpful**.

*   **B) Confocal fluorescence microscopy after Mito-RTP staining of the cells:**
    *   **Careful Point 1: Identify the reagent's purpose.** The prefix "Mito-" indicates the probe is targeted to the mitochondria. The candidates suggest "RTP" could stand for a probe measuring temperature, reactive oxygen species (ROS), or pyrophosphate.
    *   **Careful Point 2: Relate the measurement to function.** All these possibilities (temperature, ROS, pyrophosphate) are directly related to mitochondrial metabolic activity. Heat is a byproduct of respiration, ROS are generated by the electron transport chain, and pyrophosphate is involved in metabolic reactions.
    *   **Conclusion:** Using a mitochondria-specific probe to measure a parameter related to metabolic activity or stress is a valid approach. This experiment **is helpful**.

*   **C) Transformation of cells with recombinant luciferase and luminometer reading after 5 μM of luciferin addition to the supernatant:**
    *   **Careful Point 1: Understand the assay.** The luciferase-luciferin reaction produces light in a reaction that requires **ATP**. The amount of light is proportional to the ATP concentration.
    *   **Careful Point 2: Relate ATP to mitochondria.** As stated in the prompt, mitochondria are the cell's primary ATP producers. Measuring ATP levels is a fundamental way to assess the primary output of mitochondrial function.
    *   **Careful Point 3: Consider the sample location ("supernatant").** While measuring intracellular ATP is more direct, measuring extracellular ATP can be an indicator of cell stress or death, which can be triggered by mitochondrial toxins. Therefore, it is still an indirect but potentially useful readout.
    *   **Conclusion:** This experiment measures the end-product of mitochondrial energy production and **is helpful**.

*   **D) Differential centrifugation extraction of mitochondria followed by the Glucose Uptake Colorimetric Assay Kit:**
    *   **Careful Point 1: Understand the procedure.** This experiment first isolates mitochondria from the rest of the cell and then attempts to measure glucose uptake in this isolated fraction.
    *   **Careful Point 2: Analyze the biological premise.** Glucose uptake is the process of transporting glucose from outside the cell, across the plasma membrane, and into the cytoplasm. This is performed by glucose transporters (GLUTs) on the cell surface. The subsequent breakdown of glucose (glycolysis) also occurs in the cytoplasm.
    *   **Careful Point 3: Identify the flaw.** Mitochondria do **not** take up glucose. They reside inside the cell and import the products of glycolysis, primarily pyruvate. The machinery for glucose uptake is not present on mitochondrial membranes.
    *   **Conclusion:** This experimental design is based on a fundamentally incorrect biological premise. It attempts to measure a process (glucose uptake) in an organelle (the mitochondrion) that does not perform it. Therefore, this experiment **will not be helpful**.

**Final Decision**
Experiments A, B, and C all measure valid parameters of mitochondrial function or health: membrane potential, metabolic activity/stress, and ATP levels, respectively. Experiment D is fundamentally flawed because mitochondria do not perform glucose uptake.

<<<D>>>
"""

# Run the check
result = check_answer_correctness(llm_final_answer)
print(result)