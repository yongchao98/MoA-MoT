import re

def check_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given immunology question.

    The function simulates the reasoning process by:
    1. Defining the key facts presented in the question.
    2. Defining the properties of each possible answer (immunological process).
    3. Comparing the LLM's chosen answer against the facts to determine its validity.
    """

    # --- Step 1: Parse the LLM's answer ---
    # The final answer is expected in the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure to parse the answer. The answer format should be <<<A>>>, <<<B>>>, etc."
    
    llm_choice = match.group(1)

    # --- Step 2: Define the question's constraints and options ---
    # The options are mapped based on the final provided answer block.
    options_map = {
        'A': 'class switching recombination',
        'B': 'somatic hypermutation',
        'C': 'VDJ recombination',
        'D': 'complement activation'
    }
    
    llm_process_name = options_map.get(llm_choice)
    if not llm_process_name:
        return f"Invalid option '{llm_choice}'. The options are A, B, C, D."

    # Key facts derived from the question text
    question_facts = {
        "location": "secondary_lymphoid_organ",  # Peyer's patches
        "timing": "post_antigen_exposure",        # Response to antigen delivery
        "cell_state": "proliferating",            # "proliferating cell population"
        "genetic_target": "variable_region",      # "variable heavy chain gene"
        "genetic_effect": "high_variability"      # "high variability"
    }

    # --- Step 3: Define the properties of each immunological process ---
    process_properties = {
        "somatic hypermutation": {
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_exposure",
            "genetic_target": "variable_region",
            "genetic_effect": "high_variability",
            "is_correct": True
        },
        "class switching recombination": {
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_exposure",
            "genetic_target": "constant_region",  # Mismatch with question
            "genetic_effect": "isotype_change",
            "is_correct": False,
            "reason": "Class switching recombination affects the *constant* region of the heavy chain, not the *variable* region where high variability was observed."
        },
        "VDJ recombination": {
            "location": "primary_lymphoid_organ",  # Mismatch with question
            "timing": "pre_antigen_exposure",     # Mismatch with question
            "genetic_target": "variable_region",
            "genetic_effect": "initial_diversity",
            "is_correct": False,
            "reason": "VDJ recombination occurs in primary lymphoid organs (e.g., bone marrow) *before* antigen encounter, not in secondary lymphoid organs *after* antigen stimulation as described."
        },
        "complement activation": {
            "location": "plasma/tissues",
            "timing": "post_antigen_exposure",
            "genetic_target": "none",             # Mismatch with question
            "genetic_effect": "none",
            "is_correct": False,
            "reason": "Complement activation is a protein cascade and does not involve genetic changes in B cells, so it cannot explain the observed genetic variability."
        }
    }

    # --- Step 4: Evaluate the LLM's choice ---
    chosen_process_info = process_properties[llm_process_name]

    if chosen_process_info["is_correct"]:
        # The LLM chose the process that matches all facts.
        return "Correct"
    else:
        # The LLM chose an incorrect process. Return the reason why it's wrong.
        return f"Incorrect. The selected answer was '{llm_process_name}'. Reason: {chosen_process_info['reason']}"

# The user-provided question and the final answer to be checked.
question = """
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) class switching recombination
B) somatic hypermutation
C) VDJ recombination
D) complement activation
"""

llm_answer = """
Here is a step-by-step analysis of the question to determine the correct process.

1.  **Deconstruct the Experimental Scenario:**
    *   **Antigen and Delivery:** Rotavirus proteins (antigens) are delivered orally, targeting the gut's immune system.
    *   **Location of Response:** The response is observed in the Peyer's patches. These are secondary lymphoid organs in the intestine, which are primary sites for initiating immune responses to gut antigens.
    *   **Cellular Activity:** A "proliferating cell population" is isolated. This indicates that lymphocytes (in this case, B cells) have been activated by the antigen and are undergoing clonal expansion. This proliferation occurs in specialized structures called germinal centers within the Peyer's patches.
    *   **Genetic Finding:** The "variable heavy chain gene" is sequenced, and "high variability" is found. This is the most critical piece of information. The variable region of an antibody is what binds to the antigen.

2.  **Evaluate the Options Based on the Evidence:**

    *   **A) Class switching recombination (CSR):** This process occurs in activated B cells within germinal centers. However, it modifies the **constant** region of the heavy chain gene, not the variable region. This changes the antibody's function (e.g., switching from IgM to IgA for mucosal immunity) but does not alter its antigen-binding site or explain the variability observed there.

    *   **B) Somatic hypermutation (SHM):** This process occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs. Its specific function is to introduce a high rate of point mutations into the genes coding for the **variable** regions of the antibody heavy and light chains. This creates a pool of B cells with slightly different antigen receptors. Those with higher affinity for the antigen are selected to survive and proliferate further (affinity maturation). This process directly and perfectly explains the "high variability" observed in the "variable heavy chain gene" of the proliferating cell population.

    *   **C) VDJ recombination:** This is the genetic process that creates the initial, vast diversity of B cell receptors. However, it occurs in the bone marrow during the development of B cells, **before** they ever encounter an antigen. The experiment describes a process happening in a secondary lymphoid organ **after** antigen stimulation. Therefore, VDJ recombination is not the process taking place at this stage.

    *   **D) Complement activation:** This is a system of plasma proteins that is part of the innate immune response. It helps clear pathogens but does not involve any genetic changes in B cells. It is unrelated to the observed genetic variability.

3.  **Conclusion:**
    The combination of location (Peyer's patch), cellular activity (proliferation after antigen exposure), and the specific genetic outcome (high variability in the variable heavy chain gene) are the classic hallmarks of somatic hypermutation.

<<<B>>>
"""

# Run the check
result = check_correctness(question, llm_answer)
print(result)