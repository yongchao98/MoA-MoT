import re

def check_correctness_of_answer(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer for the immunology question.

    The function models the key characteristics of the immunological processes and compares them
    against the conditions described in the experiment to determine the correct answer. It then
    evaluates the LLM's provided answer.

    Args:
        llm_response: A string containing the LLM's full response, including the final answer
                      in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness of the answer.
    """

    # Define the characteristics of each immunological process based on established science.
    # The lettering (A, B, C, D) is based on the final analysis provided in the prompt.
    processes = {
        "A": {
            "name": "VDJ recombination",
            "location": "primary_lymphoid_organ",  # e.g., Bone marrow
            "timing": "pre_antigen_encounter",
            "target_region": "variable_region",
            "description": "Creates the initial diversity of B-cell receptors before antigen exposure."
        },
        "B": {
            "name": "Complement activation",
            "location": "plasma/tissues",
            "timing": "post_antigen_encounter",
            "target_region": "N/A",  # Not a genetic modification of B-cells
            "description": "A protein cascade that helps clear pathogens; does not alter B-cell genes."
        },
        "C": {
            "name": "Class switching recombination",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "target_region": "constant_region",  # Key distinction
            "description": "Changes the antibody's function (isotype) by altering the constant region."
        },
        "D": {
            "name": "Somatic hypermutation",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "target_region": "variable_region",
            "description": "Introduces point mutations into the variable region to increase antibody affinity."
        }
    }

    # Define the key conditions observed in the experiment described in the question.
    experiment_conditions = {
        "location": "secondary_lymphoid_organ",  # Peyer's patches
        "timing": "post_antigen_encounter",      # Response to the delivered rotavirus protein
        "target_region": "variable_region",     # "high variability in the variable heavy chain gene"
    }

    # Determine the correct answer by finding the process that matches all experimental conditions.
    correct_option = None
    for option, details in processes.items():
        if (details.get("location") == experiment_conditions["location"] and
            details.get("timing") == experiment_conditions["timing"] and
            details.get("target_region") == experiment_conditions["target_region"]):
            correct_option = option
            break

    # Extract the LLM's final answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Error: Could not find the final answer in the required format '<<<X>>>' in the provided text."

    llm_option = match.group(1)

    # Compare the LLM's answer with the determined correct answer.
    if llm_option == correct_option:
        return "Correct"
    else:
        # If incorrect, construct a detailed reason.
        llm_choice_details = processes.get(llm_option)
        correct_choice_details = processes.get(correct_option)

        reason = (f"The provided answer is '{llm_option}' ({llm_choice_details['name']}), which is incorrect. "
                  f"The correct answer is '{correct_option}' ({correct_choice_details['name']}).\n\n"
                  f"Reasoning:\n")

        # Explain why the LLM's choice is wrong by checking against the experimental conditions.
        if llm_choice_details.get("location") != experiment_conditions["location"]:
            reason += (f"- The chosen process, {llm_choice_details['name']}, occurs in a '{llm_choice_details.get('location')}', "
                       f"but the experiment takes place in a '{experiment_conditions['location']}' (Peyer's patches).\n")
        
        if llm_choice_details.get("timing") != experiment_conditions["timing"]:
            reason += (f"- The chosen process, {llm_choice_details['name']}, occurs '{llm_choice_details.get('timing')}', "
                       f"but the experiment observes a response '{experiment_conditions['timing']}'.\n")
        
        if llm_choice_details.get("target_region") != experiment_conditions["target_region"]:
            reason += (f"- The chosen process, {llm_choice_details['name']}, affects the '{llm_choice_details.get('target_region')}', "
                       f"but the experiment observed high variability in the '{experiment_conditions['target_region']}'.\n")
            
        reason += (f"\nThe correct process, {correct_choice_details['name']}, perfectly matches all conditions: "
                   f"it occurs '{experiment_conditions['timing']}' in a '{experiment_conditions['location']}' "
                   f"and introduces high variability into the '{experiment_conditions['target_region']}'.")
        
        return reason

# The final response from the LLM to be checked
llm_final_response = """
Here is a step-by-step analysis of the question to determine the correct process.

**Step 1: Deconstruct the experimental scenario.**
*   **Antigen and Delivery:** Rotavirus inner capsid proteins (an antigen) are delivered orally. This targets the mucosal immune system of the gastrointestinal tract.
*   **Location of Immune Response:** The response is observed in the Peyer's patches. These are secondary lymphoid organs located in the small intestine, which are primary sites for initiating immune responses to gut antigens.
*   **Cellular Activity:** A "proliferating cell population" is isolated. This indicates that lymphocytes (in this case, B cells) have been activated by the antigen and are undergoing clonal expansion. This proliferation typically occurs in specialized structures called germinal centers within the Peyer's patches.
*   **Genetic Analysis (Key Finding):** The "variable heavy chain gene" of these proliferating cells is sequenced, and "high variability" is observed. The variable region of the heavy chain is the part of the antibody that binds to the antigen.

**Step 2: Evaluate the options based on the key finding.**
The central task is to identify the process that generates high variability in the *variable region* of immunoglobulin genes in *proliferating B cells* within a *secondary lymphoid organ* *after* antigen exposure.

*   **A) VDJ recombination:** This process creates the initial diversity of B cell receptors. However, it occurs in the bone marrow (a primary lymphoid organ) during the development of B cells, *before* they encounter an antigen. The scenario describes a process happening after antigen encounter in a secondary lymphoid organ. Therefore, VDJ recombination is incorrect.

*   **B) Complement activation:** This is a system of plasma proteins that helps clear pathogens. It is a component of the innate and adaptive immune response but does not involve any genetic changes to the B cells themselves. Therefore, complement activation is incorrect.

*   **C) Class switching recombination (CSR):** This process also occurs in activated B cells within germinal centers. However, CSR modifies the *constant* region of the heavy chain gene, not the *variable* region. This changes the antibody's function (e.g., switching from IgM to IgA) but does not alter the antigen-binding site or create variability there. Therefore, class switching recombination is not the direct cause of the observed variability.

*   **D) Somatic hypermutation (SHM):** This is a process that occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs (like Peyer's patches). SHM specifically introduces a high rate of point mutations into the genes encoding the *variable regions* of the heavy and light chains. This creates a pool of B cells with slightly different antibodies, allowing for the selection of those with higher affinity for the antigen (affinity maturation). This process perfectly explains the "high variability" observed in the "variable heavy chain gene" of the "proliferating cell population."

**Step 3: Conclude the analysis.**
The experimental evidence—proliferating B cells in a Peyer's patch showing high genetic variability specifically in the variable region of their immunoglobulin genes after antigen stimulation—is the classic signature of somatic hypermutation.

<<<D>>>
"""

# Run the check
result = check_correctness_of_answer(llm_final_response)
print(result)