import re

def check_immunology_answer(llm_answer_text):
    """
    Checks the correctness of the final answer for the immunology question.

    The function verifies the answer based on established immunological principles
    related to B cell development and activation.
    """

    # --- Define Ground Truth ---

    # 1. Key constraints derived from the question's experimental setup
    question_constraints = {
        "location": "secondary lymphoid organ",  # Peyer's patch
        "timing": "post-antigen stimulation",   # Proliferating cells after antigen delivery
        "genetic_target": "variable region",    # "variable heavy chain gene"
        "observation": "high variability"       # The key finding
    }

    # 2. Properties of the possible processes (the options)
    processes = {
        "class switching recombination": {
            "location": "secondary lymphoid organ",
            "timing": "post-antigen stimulation",
            "genetic_target": "constant region",  # Mismatch with question
            "observation": "change in antibody isotype"
        },
        "somatic hypermutation": {
            "location": "secondary lymphoid organ",
            "timing": "post-antigen stimulation",
            "genetic_target": "variable region",    # Match
            "observation": "high variability"       # Match
        },
        "VDJ recombination": {
            "location": "primary lymphoid organ",   # Mismatch with question
            "timing": "pre-antigen stimulation",    # Mismatch with question
            "genetic_target": "variable region",
            "observation": "creation of initial receptor diversity"
        },
        "complement activation": {
            "location": "blood/tissues",
            "timing": "post-antigen stimulation",
            "genetic_target": "N/A (not a B-cell genetic process)", # Mismatch
            "observation": "protein cascade for pathogen clearance"
        }
    }
    
    # 3. The options as presented in the original question
    original_options = {
        "A": "class switching recombination",
        "B": "somatic hypermutation",
        "C": "VDJ recombination",
        "D": "complement activation"
    }

    # --- Logic to Determine Correct Answer ---

    correct_process_name = None
    for process, properties in processes.items():
        if (properties["location"] == question_constraints["location"] and
            properties["timing"] == question_constraints["timing"] and
            properties["genetic_target"] == question_constraints["genetic_target"] and
            properties["observation"] == question_constraints["observation"]):
            correct_process_name = process
            break
    
    if not correct_process_name:
        return "Error in checker logic: Could not determine the correct process."

    correct_option_letter = None
    for letter, text in original_options.items():
        if text == correct_process_name:
            correct_option_letter = letter
            break

    if not correct_option_letter:
        return f"Error in checker logic: The correct process '{correct_process_name}' was not found in the options."

    # --- Evaluate the LLM's Answer ---

    # Extract the final letter choice from the provided answer text
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect format: The answer does not contain the '<<<X>>>' format at the end."

    llm_choice_letter = match.group(1)

    # Compare and generate feedback
    if llm_choice_letter == correct_option_letter:
        return "Correct"
    else:
        llm_chosen_process = original_options.get(llm_choice_letter, "an invalid option")
        reason = (
            f"Incorrect. The provided answer is '{llm_choice_letter}', which corresponds to '{llm_chosen_process}'.\n"
            f"The correct answer is '{correct_option_letter}', which is '{correct_process_name}'.\n\n"
            "Reasoning:\n"
            "1.  **Location and Timing:** The process occurs in a Peyer's patch (a secondary lymphoid organ) in a 'proliferating cell population' after antigen stimulation. This rules out VDJ recombination, which happens in a primary lymphoid organ (bone marrow) before antigen encounter.\n"
            "2.  **Genetic Target:** The key finding is 'high variability' in the 'variable heavy chain gene'.\n"
            "3.  **Conclusion:** Somatic hypermutation is the only process that specifically introduces high rates of mutation into the *variable* region of immunoglobulin genes in proliferating B cells within a secondary lymphoid organ. Class switching recombination, while occurring at the same time and place, affects the *constant* region, not the variable region."
        )
        return reason

# The final answer from the prompt to be checked
llm_final_answer = """
Here is a step-by-step analysis of the question, incorporating careful points to distinguish between the possible answers.

**Step 1: Deconstruct the Experimental Scenario**
The experiment provides several key pieces of information that point to a specific immunological process:
*   **Antigen and Delivery:** Rotavirus proteins are delivered orally. This specifically targets the Gut-Associated Lymphoid Tissue (GALT), which is a major part of the mucosal immune system.
*   **Location of Response:** The response occurs in the Peyer's patches. These are well-defined secondary lymphoid organs in the small intestine, which are primary sites for initiating adaptive immune responses to gut antigens.
*   **Cellular Response:** A "proliferating cell population" is isolated. This indicates that lymphocytes that recognize the antigen have been activated and are undergoing clonal expansion. This rapid proliferation of B cells occurs in specialized microenvironments called germinal centers, which form within secondary lymphoid organs like Peyer's patches.
*   **Genetic Analysis and Key Finding:** The "variable heavy chain gene" is sequenced, and "high variability" is found. This is the most critical clue. It tells us:
    1.  The proliferating cells are B cells (as they possess heavy chain genes).
    2.  The genetic changes are happening in the *variable* region, which is the part of the antibody that binds to the antigen.
    3.  "High variability" means that many new point mutations have been introduced into this gene *during* the proliferation phase.

**Step 2: Evaluate the Options with Careful Distinctions**

The question asks to identify the process responsible for creating high variability in the variable region of the heavy chain gene *after* antigen stimulation in a secondary lymphoid organ.

*   **A) class switching recombination (CSR)**
    *   **Careful Point:** While CSR does occur in activated, proliferating B cells in germinal centers (correct timing and location), it exclusively modifies the **constant region** of the heavy chain gene. This changes the antibody's function (e.g., from IgM to IgA for a mucosal response) but does **not** alter the variable region or the antigen-binding site. Therefore, it cannot explain the observed variability in the variable region.

*   **B) somatic hypermutation (SHM)**
    *   **Careful Point:** This process perfectly matches all experimental observations. It occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs (like Peyer's patches). Its specific function is to introduce a high rate of point mutations into the **variable regions** of the immunoglobulin heavy and light chain genes. This is the direct cause of the "high variability" observed and is the mechanism for affinity maturation, where B cells producing higher-affinity antibodies are selected.

*   **C) VDJ recombination**
    *   **Careful Point:** This is the process that creates the initial, primary diversity of B cell receptors. The crucial distinction is that VDJ recombination occurs during B cell development in a **primary lymphoid organ** (the bone marrow) **before** the B cell ever encounters an antigen. The experiment is observing a process happening in a **secondary lymphoid organ** **after** antigen encounter.

*   **D) complement activation**
    *   **Careful Point:** Complement is a system of proteins in the blood and tissues that is part of the humoral immune response. It helps clear pathogens but does not involve any genetic modification of B cells. It cannot be the cause of genetic variability in a B cell population.

**Step 3: Conclusion**
The combination of location (Peyer's patch), timing (post-antigen stimulation), cellular activity (proliferation), and specific genetic outcome (high variability in the variable heavy chain gene) is the definitive signature of **somatic hypermutation**.

<<<B>>>
"""

# Run the check
result = check_immunology_answer(llm_final_answer)
print(result)