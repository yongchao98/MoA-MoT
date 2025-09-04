import re

def check_immunology_answer(llm_response: str) -> str:
    """
    Checks the correctness of the answer to the immunology question.

    The function verifies the answer by encoding the key facts of the experimental
    scenario and the defining features of each immunological process. It then
    logically determines the correct answer and compares it to the provided one.

    Args:
        llm_response: The full text of the LLM's response, which should
                      contain the final answer in the format <<<X>>>.

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """
    # 1. Define the key features from the experimental scenario in the question
    scenario_features = {
        "location": "secondary lymphoid organ",  # Peyer's patches
        "timing": "post-antigen encounter",      # Response to antigen delivery
        "cellular_activity": "proliferation",    # Proliferating cell population
        "genetic_target": "variable region",     # "variable heavy chain gene"
        "genetic_effect": "high variability"     # "high variability" from point mutations
    }

    # 2. Define the properties of the immunological processes given as options
    processes = {
        "A": {
            "name": "complement activation",
            "description": "A protein-based cascade in the innate immune system; does not involve genetic changes in B cells."
        },
        "B": {
            "name": "VDJ recombination",
            "timing": "pre-antigen encounter",
            "location": "primary lymphoid organ", # Bone marrow
            "genetic_target": "variable region",
            "description": "Creates initial receptor diversity before antigen encounter."
        },
        "C": {
            "name": "somatic hypermutation",
            "timing": "post-antigen encounter",
            "location": "secondary lymphoid organ",
            "cellular_activity": "proliferation",
            "genetic_target": "variable region",
            "genetic_effect": "high variability",
            "description": "Introduces point mutations in the variable region after antigen encounter to refine affinity."
        },
        "D": {
            "name": "class switching recombination",
            "timing": "post-antigen encounter",
            "location": "secondary lymphoid organ",
            "cellular_activity": "proliferation",
            "genetic_target": "constant region", # This is the key difference from SHM
            "description": "Changes the antibody's function by altering the constant region, not its antigen-binding site."
        }
    }

    # 3. Determine the correct option by matching scenario features
    correct_option = None
    for option, properties in processes.items():
        # Check if all features from the scenario match the process properties
        is_match = all(
            properties.get(key) == value for key, value in scenario_features.items()
        )
        if is_match:
            correct_option = option
            break
    
    if not correct_option:
        return "Error in checker logic: Could not programmatically determine the correct answer based on the defined features."

    # 4. Extract the LLM's answer from the response text
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Incorrect format: The final answer in the format <<<X>>> was not found in the response."

    llm_choice = match.group(1)

    # 5. Compare the LLM's choice with the derived correct answer and return the result
    if llm_choice == correct_option:
        return "Correct"
    else:
        incorrect_process = processes[llm_choice]
        correct_process = processes[correct_option]
        
        # Generate a specific reason for the error
        reason = f"The provided answer '{llm_choice}' ({incorrect_process['name']}) is incorrect. "
        
        if incorrect_process.get("genetic_target") != scenario_features["genetic_target"]:
            reason += f"The question specifies high variability in the *{scenario_features['genetic_target']}*. {incorrect_process['name']} affects the *{incorrect_process.get('genetic_target')}*."
        elif incorrect_process.get("timing") != scenario_features["timing"]:
            reason += f"The question describes a process happening *{scenario_features['timing']}*. {incorrect_process['name']} occurs *{incorrect_process.get('timing')}*."
        else:
            reason += f"The correct answer is '{correct_option}' ({correct_process['name']}) because it is the only process that matches all conditions: occurring post-antigen encounter in proliferating B cells in a secondary lymphoid organ to cause high variability in the variable gene region."
            
        return reason

# The user-provided response to be checked
llm_answer_text = """
Here is a plan to solve the problem.
*   Analyze the experimental scenario described in the question to identify the key immunological events.
*   Evaluate each of the four options (complement activation, VDJ recombination, somatic hypermutation, class switching recombination) against the experimental evidence.
*   Review the provided candidate answers to find the consensus in their reasoning, paying close attention to how they justify their choice.
*   Synthesize the analysis to determine the single best answer and map it to the correct letter from the original question's options.

**Step 1: Analyze the Experimental Scenario**
The question describes a classic immunology experiment. Let's break down the key points:
*   **Antigen Delivery:** Oral gavage of rotavirus proteins targets the Gut-Associated Lymphoid Tissue (GALT).
*   **Location of Response:** The response is observed in Peyer's patches, which are secondary lymphoid organs in the intestine designed to respond to gut antigens.
*   **Cellular Activity:** A "proliferating cell population" is found. This indicates clonal expansion of lymphocytes that have recognized the antigen. In secondary lymphoid organs, this occurs in specialized structures called germinal centers.
*   **Genetic Finding:** Sequencing reveals "high variability" specifically in the "variable heavy chain gene." This gene region codes for the antigen-binding site of antibodies, and the term "heavy chain" confirms the proliferating cells are B cells.

The core finding is the generation of high genetic diversity in the antigen-binding region of antibodies in a population of proliferating B cells, *after* they have been stimulated by an antigen in a secondary lymphoid organ.

**Step 2: Evaluate the Options**
Let's assess which process matches these findings:
*   **A) complement activation:** This is a system of proteins in the blood that helps clear pathogens. It is an effector mechanism of the immune system and does not involve any genetic changes in B cells. This is incorrect.
*   **B) VDJ recombination:** This is the genetic process that creates the initial diversity of B cell receptors. However, it occurs in the bone marrow (a primary lymphoid organ) during B cell development, *before* the cell ever encounters an antigen. The experiment describes a process happening *after* antigen encounter in a secondary lymphoid organ. This is incorrect.
*   **C) somatic hypermutation (SHM):** This process occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs (like Peyer's patches). Its specific function is to introduce a high rate of point mutations into the genes coding for the *variable regions* of antibodies. This creates the "high variability" observed in the experiment and allows for the selection of B cells with higher affinity for the antigen. This perfectly matches all aspects of the scenario.
*   **D) class switching recombination (CSR):** This process also occurs in activated B cells in germinal centers. However, it affects the *constant* region of the heavy chain, not the *variable* region. It changes the antibody's function (e.g., switching from IgM to IgA for mucosal immunity) but does not alter the antigen-binding site or create the variability observed there. This is incorrect.

**Step 3: Review Candidate Answers**
All 17 candidate answers, despite some having different final letter choices due to variations in their option ordering, present a unanimous line of reasoning. They all correctly identify that:
*   The location (Peyer's patch), cell activity (proliferation), and timing (post-antigen exposure) point to a germinal center reaction.
*   The key finding of "high variability" in the "variable heavy chain gene" is the defining characteristic of somatic hypermutation.
*   They all correctly rule out the other options for the reasons outlined in Step 2.

**Step 4: Synthesize and Conclude**
The experimental evidence points unequivocally to somatic hypermutation. The scenario describes the precise location, cellular activity, and genetic outcome of this process. It is the mechanism by which the adaptive immune system refines the quality of the antibody response after an initial antigen encounter.

<<<C>>>
"""

# Run the check
result = check_immunology_answer(llm_answer_text)
print(result)