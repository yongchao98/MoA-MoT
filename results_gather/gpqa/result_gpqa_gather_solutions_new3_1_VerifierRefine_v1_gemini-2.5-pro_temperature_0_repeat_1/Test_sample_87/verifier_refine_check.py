import re

def check_answer(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer to the immunology question.

    The function formalizes the key constraints from the question and checks if the
    selected answer satisfies all of them.
    """

    # Define the key facts/constraints from the question
    question_constraints = {
        "location": "secondary lymphoid organ",  # Peyer's patches
        "timing": "post-antigen encounter",     # After oral gavage and uptake
        "cell_state": "proliferating",          # Proliferating cell population
        "genetic_target": "variable region",    # "variable heavy chain gene"
        "genetic_effect": "high variability"    # "high variability" observed
    }

    # Define the properties of each immunological process
    processes = {
        "complement activation": {
            "location": None,
            "timing": None,
            "cell_state": None,
            "genetic_target": None,
            "genetic_effect": "protein cascade",
            "is_correct": False,
            "reason": "It is a protein cascade, not a process of genetic mutation in B cells."
        },
        "class switching recombination": {
            "location": "secondary lymphoid organ",
            "timing": "post-antigen encounter",
            "cell_state": "proliferating",
            "genetic_target": "constant region", # This is the key difference
            "genetic_effect": "isotype switch",
            "is_correct": False,
            "reason": "It affects the constant region of the heavy chain, not the variable region as stated in the question."
        },
        "VDJ recombination": {
            "location": "primary lymphoid organ", # Key difference
            "timing": "pre-antigen encounter",   # Key difference
            "cell_state": "developing",
            "genetic_target": "variable region",
            "genetic_effect": "initial diversity generation",
            "is_correct": False,
            "reason": "It occurs in primary lymphoid organs (bone marrow) before antigen encounter, not in secondary organs after antigen stimulation."
        },
        "somatic hypermutation": {
            "location": "secondary lymphoid organ",
            "timing": "post-antigen encounter",
            "cell_state": "proliferating",
            "genetic_target": "variable region",
            "genetic_effect": "high variability",
            "is_correct": True,
            "reason": "This process perfectly matches all conditions: it occurs in proliferating B cells in secondary lymphoid organs after antigen encounter and causes high variability in the variable region."
        }
    }

    # Map the options from the question text to the process names
    # A) complement activation, B) class switching recombination, C) VDJ recombination, D) somatic hypermutation
    option_mapping = {
        'A': 'complement activation',
        'B': 'class switching recombination',
        'C': 'VDJ recombination',
        'D': 'somatic hypermutation'
    }

    # Extract the final letter answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."

    llm_choice_letter = match.group(1)
    llm_chosen_process_name = option_mapping.get(llm_choice_letter)

    if not llm_chosen_process_name:
        return f"Incorrect: The chosen letter '{llm_choice_letter}' does not correspond to any option."

    # Check if the chosen process is the correct one
    chosen_process_details = processes[llm_chosen_process_name]

    if chosen_process_details["is_correct"]:
        # Further check if the LLM's reasoning aligns with the facts
        reasoning_text = llm_answer_text.lower()
        if "somatic hypermutation" not in reasoning_text and "shm" not in reasoning_text:
             return "Incorrect: The final letter choice is correct, but the reasoning does not mention the correct process (somatic hypermutation)."
        if "variable region" not in reasoning_text:
             return "Incorrect: The final letter choice is correct, but the reasoning fails to mention the key detail that the 'variable region' is affected."
        if "secondary lymphoid organ" not in reasoning_text and "peyer's patch" not in reasoning_text:
             return "Incorrect: The final letter choice is correct, but the reasoning fails to mention the correct location (secondary lymphoid organ / Peyer's patch)."
        if "post-antigen" not in reasoning_text and "after antigen" not in reasoning_text:
             return "Incorrect: The final letter choice is correct, but the reasoning fails to mention the correct timing (post-antigen encounter)."
        
        return "Correct"
    else:
        # The LLM chose the wrong answer, provide the reason why it's wrong
        correct_process_letter = [k for k, v in option_mapping.items() if processes[v]["is_correct"]][0]
        return (f"Incorrect: The chosen answer '{llm_choice_letter}' ({llm_chosen_process_name}) is wrong. "
                f"Reason: {chosen_process_details['reason']} "
                f"The correct answer is '{correct_process_letter}' (somatic hypermutation), which is the only process that causes high variability in the variable region of proliferating B cells in a secondary lymphoid organ after antigen encounter.")


# The user-provided question and answer to check
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) complement activation
B) class switching recombination
C) VDJ recombination
D) somatic hypermutation
"""

llm_answer = """
Here is a step-by-step analysis of the question and the provided options.

**Step 1: Deconstruct the Experimental Scenario**
*   **Antigen & Delivery:** GFP-tagged rotavirus proteins are delivered orally. This targets the mucosal immune system, specifically the Gut-Associated Lymphoid Tissue (GALT).
*   **Location of Response:** The response is observed in Peyer's patches, which are secondary lymphoid organs in the small intestine. They are primary sites for initiating immune responses to antigens from the gut lumen. The fluorescence confirms the antigen has reached these sites.
*   **Cellular Response:** A "proliferating cell population" is isolated. This indicates clonal expansion, a key feature of an adaptive immune response where lymphocytes that recognize a specific antigen multiply rapidly.
*   **Genetic Analysis:** The "variable heavy chain gene" of these proliferating cells is sequenced. This points specifically to B cells, as this gene codes for the antigen-binding portion of the B cell receptor (and secreted antibodies).
*   **Key Finding:** The sequencing reveals "high variability" in this gene. This means that among the population of proliferating B cells, there are numerous different mutations in the gene segment that codes for the antigen-binding site.

**Step 2: Evaluate the Immunological Processes (The Options)**

The question asks to identify the process responsible for creating this high variability in the variable region of the heavy chain gene *after* antigen stimulation in a secondary lymphoid organ.

*   **A) Complement activation:** This is a system of plasma proteins that is part of the innate and humoral immune response. It helps clear pathogens but does not involve any genetic changes in B cells. Therefore, it cannot be the cause of the observed genetic variability.

*   **B) Class switching recombination (CSR):** This process occurs in activated B cells within germinal centers (found in Peyer's patches). However, CSR modifies the *constant* region of the heavy chain gene, not the *variable* region. This changes the antibody isotype (e.g., from IgM to IgA), which alters its function, but it does not affect the antigen-binding site or its variability.

*   **C) V(D)J recombination:** This is the process that creates the initial, vast diversity of B cell receptors. Crucially, this occurs during the development of B cells in the bone marrow, *before* they ever encounter an antigen. The experiment is observing a process happening in mature, proliferating B cells in a secondary lymphoid organ *after* antigen encounter. Therefore, V(D)J recombination is not the process taking place at this stage.

*   **D) Somatic hypermutation (SHM):** This process occurs in activated B cells that are proliferating within germinal centers of secondary lymphoid organs (like Peyer's patches). SHM specifically introduces a high rate of point mutations into the genes encoding the *variable* regions of the antibody heavy and light chains. This creates a pool of B cells with slightly different antigen receptors. B cells whose receptors have a higher affinity for the antigen are then selected to survive and proliferate. This process, known as affinity maturation, directly and perfectly explains the "high variability" observed in the "variable heavy chain gene" of the proliferating B cell population.

**Step 3: Conclusion**

The combination of the location (Peyer's patch), the cell type (proliferating B cells), the timing (after antigen stimulation), and the specific genetic outcome (high variability in the variable region) are all classic hallmarks of somatic hypermutation.

<<<D>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)