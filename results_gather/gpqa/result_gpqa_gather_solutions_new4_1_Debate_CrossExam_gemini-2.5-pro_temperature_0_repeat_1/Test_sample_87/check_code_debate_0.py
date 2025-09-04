import re

def check_immunology_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given immunology question.

    The function verifies if the chosen option correctly identifies the process based on
    the key constraints provided in the question:
    1. Location: Secondary lymphoid organ (Peyer's patch)
    2. Timing: Post-antigen encounter (in response to the virus protein)
    3. Cellular Activity: Proliferation
    4. Genetic Target: Variable region of the immunoglobulin gene
    5. Genetic Outcome: High variability (point mutations)
    """

    # Define the properties of each immunological process based on the canonical options
    # A) class switching recombination, B) complement activation, C) VDJ recombination, D) somatic hypermutation
    processes = {
        "A": {
            "name": "Class switching recombination",
            "genetic_target": "constant_region",  # CRITICAL MISMATCH: Affects constant, not variable region.
            "timing": "post_antigen_encounter",
            "location": "secondary_lymphoid_organ"
        },
        "B": {
            "name": "Complement activation",
            "type": "protein_cascade",  # CRITICAL MISMATCH: Not a genetic modification process.
            "genetic_target": None
        },
        "C": {
            "name": "VDJ recombination",
            "genetic_target": "variable_region",
            "timing": "pre_antigen_encounter",  # CRITICAL MISMATCH: Occurs before antigen encounter.
            "location": "primary_lymphoid_organ"  # CRITICAL MISMATCH: Occurs in primary, not secondary organs.
        },
        "D": {
            "name": "Somatic hypermutation",
            "genetic_target": "variable_region",  # Match
            "timing": "post_antigen_encounter",  # Match
            "location": "secondary_lymphoid_organ",  # Match
            "cellular_activity": "proliferation",  # Match
            "genetic_outcome": "high_variability" # Match
        }
    }

    # The correct answer based on biological facts
    correct_choice = 'D'
    correct_process = processes[correct_choice]

    # Extract the final choice from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The answer format is invalid. It must end with '<<<X>>>' where X is A, B, C, or D."

    llm_choice = match.group(1)

    if llm_choice == correct_choice:
        # The LLM chose the correct letter. Now, let's quickly check if its reasoning is sound.
        # A good answer should mention the variable region, the location (Peyer's patch/secondary organ),
        # and the timing (post-antigen).
        reasoning = llm_answer_text.lower()
        if "variable region" in reasoning and "somatic hypermutation" in reasoning and ("peyer's patch" in reasoning or "secondary lymphoid" in reasoning):
             return "Correct"
        else:
             # This case is unlikely if the final letter is correct, but it's a sanity check.
             return "Correct. The final choice 'D' is correct, although the reasoning might be incomplete or missing key terms."

    else:
        # The LLM chose the wrong letter. Provide a specific reason.
        chosen_process = processes.get(llm_choice)
        if not chosen_process:
            return f"Invalid choice '{llm_choice}' provided in the answer."

        if llm_choice == 'A':
            return f"Incorrect. The chosen answer is '{chosen_process['name']}'. This process affects the *constant* region of the heavy chain, but the question specifies high variability in the *variable* region."
        elif llm_choice == 'B':
            return f"Incorrect. The chosen answer is '{chosen_process['name']}'. This is a protein-based system and does not cause genetic mutations in B cells."
        elif llm_choice == 'C':
            return f"Incorrect. The chosen answer is '{chosen_process['name']}'. This process occurs in *primary* lymphoid organs (e.g., bone marrow) *before* antigen encounter. The experiment takes place in a *secondary* lymphoid organ (Peyer's patch) *after* antigen encounter."
        else:
            return f"Incorrect. The chosen answer '{llm_choice}' does not fit the experimental evidence."


# The provided question text
question_text = """
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) class switching recombination
B) complement activation
C) VDJ recombination
D) somatic hypermutation
"""

# The LLM's final answer to be checked
llm_final_answer = """
Here is a step-by-step analysis of the question to determine the correct process.

**1. Deconstruct the Experimental Scenario**

*   **Antigen and Delivery:** GFP-tagged rotavirus proteins are delivered orally. This mimics a natural infection route and targets the gut-associated lymphoid tissue (GALT).
*   **Location of Immune Response:** The response is observed in the Peyer's patches. These are secondary lymphoid organs located in the intestine, which are primary sites for initiating immune responses to antigens from the gut.
*   **Cellular Activity:** A "proliferating cell population" is isolated. This indicates that lymphocytes (in this case, B cells) have been activated by the antigen and are undergoing clonal expansion. This proliferation occurs in specialized structures called germinal centers, which form within secondary lymphoid organs.
*   **Key Genetic Finding:** Sequencing reveals "high variability" specifically in the "variable heavy chain gene." The heavy chain is a component of the B cell receptor and antibodies. The variable region is the part that binds to the antigen. The high variability means many point mutations have been introduced into this gene among the proliferating cells.

**2. Evaluate the Options Based on the Experimental Evidence**

The core task is to identify the process that causes high variability in the variable region of immunoglobulin genes in proliferating B cells *after* they have been activated by an antigen in a secondary lymphoid organ.

*   **A) Class switching recombination:** This process occurs in activated B cells within germinal centers. However, it affects the **constant region** of the heavy chain, not the variable region. It switches the antibody's function (e.g., from IgM to IgA) but does not alter the part that binds the antigen. Therefore, it does not explain the variability in the *variable* gene.

*   **B) Complement activation:** This is a system of proteins in the blood and tissues that helps clear pathogens. It is a component of the innate immune response and is not a process that involves genetic changes in B cells. Therefore, it is incorrect.

*   **C) VDJ recombination:** This is the crucial process that creates the initial diversity of B cell receptors. However, it occurs during the development of B cells in the **primary lymphoid organs** (bone marrow), **before** the cells ever encounter an antigen. The scenario describes a process happening in a **secondary lymphoid organ** *after* antigen exposure. Therefore, this is not the correct process for this stage of the immune response.

*   **D) Somatic hypermutation (SHM):** This process is the hallmark of the germinal center reaction in secondary lymphoid organs. It is triggered after B cells are activated by an antigen. During proliferation, SHM introduces a high rate of point mutations specifically into the **variable regions** of the immunoglobulin genes. This creates a pool of B cells with slightly different receptors, allowing for the selection of those with the highest affinity for the antigen (a process called affinity maturation). This mechanism perfectly explains the observation of "high variability" in the "variable heavy chain gene" of the "proliferating cell population" in the "Peyer's patches."

**Conclusion**

The combination of location (Peyer's patch), cellular activity (proliferation after antigen encounter), and the specific genetic change (high variability in the variable region) is uniquely characteristic of somatic hypermutation.

<<<D>>>
"""

# Run the check
result = check_immunology_answer(question_text, llm_final_answer)
print(result)