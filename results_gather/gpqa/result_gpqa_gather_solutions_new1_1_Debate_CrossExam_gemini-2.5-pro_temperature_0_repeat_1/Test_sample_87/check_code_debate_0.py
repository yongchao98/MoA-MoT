import re

def check_immunology_question(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given immunology question.

    The function models the key constraints from the question and the properties of the
    four possible biological processes. It then determines the correct answer based on
    these models and compares it to the LLM's provided answer.

    Args:
        llm_answer_text: The full text of the LLM's response, which should
                         contain the final answer in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """
    # --- 1. Define Constraints from the Question ---
    # These are the facts established by the experimental scenario.
    question_constraints = {
        "location": "secondary_lymphoid_organ",  # Event occurs in Peyer's patches.
        "timing": "post_antigen_encounter",      # Occurs after antigen stimulation (proliferating cells).
        "genetic_target": "variable_region",     # High variability found in the "variable heavy chain gene".
        "is_genetic_modification": True          # Sequencing reveals genetic changes.
    }

    # --- 2. Define Properties of Each Answer Option ---
    # This dictionary models the biological facts for each process.
    processes = {
        "A": {
            "name": "complement activation",
            "is_genetic_modification": False, # It's a protein cascade, not a genetic change in lymphocytes.
            "location": "blood/tissues",
            "timing": "post_antigen_encounter",
            "genetic_target": None,
        },
        "B": {
            "name": "class switching recombination",
            "is_genetic_modification": True,
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "genetic_target": "constant_region", # This is the key mismatch.
        },
        "C": {
            "name": "somatic hypermutation",
            "is_genetic_modification": True,
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_encounter",
            "genetic_target": "variable_region", # This is a perfect match.
        },
        "D": {
            "name": "VDJ recombination",
            "is_genetic_modification": True,
            "location": "primary_lymphoid_organ", # Mismatch: occurs in bone marrow.
            "timing": "pre_antigen_encounter",   # Mismatch: occurs during B-cell development.
            "genetic_target": "variable_region",
        }
    }

    # --- 3. Determine the Logically Correct Answer ---
    correct_option = None
    for option, properties in processes.items():
        match = True
        for constraint, value in question_constraints.items():
            if properties.get(constraint) != value:
                match = False
                break
        if match:
            correct_option = option
            break

    # --- 4. Parse and Evaluate the LLM's Answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find the answer in the required format '<<<X>>>' in the provided text."

    llm_option = match.group(1)

    if llm_option == correct_option:
        return "Correct"
    else:
        # --- 5. Generate a Detailed Error Report ---
        llm_choice_props = processes[llm_option]
        correct_choice_props = processes[correct_option]
        
        reason = (f"Incorrect. The provided answer is '{llm_option}' ({llm_choice_props['name']}), "
                  f"but the correct answer is '{correct_option}' ({correct_choice_props['name']}).\n\n")
        
        # Check each constraint to explain why the LLM's choice is wrong.
        if llm_choice_props['is_genetic_modification'] != question_constraints['is_genetic_modification']:
            reason += (f"Reason: The question implies a genetic modification process to explain sequence variability. "
                       f"'{llm_choice_props['name']}' is not a genetic modification of lymphocytes.\n")
        
        if llm_choice_props['location'] != question_constraints['location']:
            reason += (f"Reason: The process '{llm_choice_props['name']}' occurs in a {llm_choice_props['location']}, "
                       f"but the question specifies the event happens in a {question_constraints['location']} (Peyer's patches).\n")

        if llm_choice_props['timing'] != question_constraints['timing']:
            reason += (f"Reason: The process '{llm_choice_props['name']}' occurs {llm_choice_props['timing']}, "
                       f"but the question describes a response that is {question_constraints['timing']}.\n")

        if llm_choice_props['genetic_target'] != question_constraints['genetic_target']:
            reason += (f"Reason: The process '{llm_choice_props['name']}' affects the {llm_choice_props['genetic_target']}, "
                       f"but the question specifies high variability in the {question_constraints['genetic_target']}.\n")

        reason += f"\nThe correct process, '{correct_choice_props['name']}', perfectly matches all conditions."
        return reason

# The provided answer from the LLM
llm_answer = """
Here is a step-by-step analysis of the question to determine the correct process.

**1. Deconstruction of the Experimental Scenario**

*   **Antigen and Delivery:** GFP-tagged rotavirus inner capsid proteins are delivered via oral gavage. This mimics a natural infection of the gut and targets the mucosal immune system.
*   **Location of Immune Response:** The antigen is found in the Peyer's patches. These are secondary lymphoid organs embedded in the small intestine, serving as the primary sites for initiating immune responses to intestinal antigens.
*   **Cellular Activity:** A "proliferating cell population" is isolated from the Peyer's patches. This indicates that lymphocytes (B or T cells) have recognized the antigen and are undergoing clonal expansion, a hallmark of an adaptive immune response.
*   **Genetic Analysis:** The "variable heavy chain gene" is sequenced. The terms "heavy chain" and "variable region" are specific to B cell receptors (BCRs) and the antibodies they produce. This strongly implies the proliferating cells are B cells.
*   **Key Observation:** The sequencing reveals "high variability" in this specific gene region among the different cells in the proliferating population.

**2. Evaluation of the Options Based on the Observations**

The core of the question is to identify the process that generates high genetic variability in the *variable region* of B cell receptor genes *after* antigen stimulation in a *secondary lymphoid organ*.

*   **A) Complement activation:** This is a system of proteins in the blood and tissues that helps clear pathogens. It is a component of the humoral immune system but does not involve any genetic modification of B cells. Therefore, it cannot be the cause of the observed genetic variability.

*   **B) Class switching recombination (CSR):** This process occurs in activated B cells within germinal centers (found in Peyer's patches). However, CSR modifies the **constant region** of the heavy chain gene, not the variable region. This changes the antibody's isotype (e.g., from IgM to IgA), which alters its function, but it does not change the antigen-binding site or create the variability observed in the *variable* region.

*   **C) Somatic hypermutation (SHM):** This is a perfect match for the observations.
    *   **Location:** It occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs, such as Peyer's patches.
    *   **Timing:** It is initiated *after* B cells are activated by an antigen.
    *   **Mechanism:** It introduces a high rate of point mutations specifically into the **variable regions** of the immunoglobulin heavy and light chain genes.
    *   **Result:** This creates the "high variability" described in the question, allowing for the selection of B cells that produce antibodies with higher affinity for the antigen (a process called affinity maturation).

*   **D) VDJ recombination:** This is the process that creates the initial, vast diversity of B cell receptors. However, it occurs in the **primary lymphoid organ** (bone marrow) during the development of B cells, **before** they ever encounter an antigen. The experiment is observing a process happening in mature, proliferating B cells in a secondary lymphoid organ *after* antigen encounter.

**3. Conclusion and Careful Points**

*   The location (Peyer's patch), timing (post-antigen stimulation), cell type (proliferating B cells), and specific genetic outcome (high variability in the variable heavy chain gene) are all classic hallmarks of **somatic hypermutation**.
*   It is crucial to distinguish SHM from VDJ recombination based on timing (response vs. development) and location (secondary vs. primary lymphoid organs).
*   It is also crucial to distinguish SHM from class switching recombination based on the genetic target (variable region vs. constant region).

The process described is unequivocally somatic hypermutation.

<<<C>>>
"""

# Run the check
result = check_immunology_question(llm_answer)
print(result)