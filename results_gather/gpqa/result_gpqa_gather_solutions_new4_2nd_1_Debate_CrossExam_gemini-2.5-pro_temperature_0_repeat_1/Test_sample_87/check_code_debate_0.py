import re

def check_answer(question, candidate_answers, final_answer_text):
    """
    Checks the correctness of the final answer for the given immunology question.

    Args:
        question (str): The question text.
        candidate_answers (list): A list of candidate answers from other models.
        final_answer_text (str): The text containing the final proposed answer.

    Returns:
        str: "Correct" if the answer is correct, or a reason for why it's incorrect.
    """

    # Extract the final answer letter from the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Invalid answer format. The final answer should be in the format <<<X>>>."
    
    proposed_answer_letter = match.group(1)

    # Define the key facts from the question scenario
    question_facts = {
        'location': 'secondary_lymphoid_organ',  # Peyer's patches
        'timing': 'post_antigen_exposure',      # Response to antigen delivery
        'cellular_process': 'proliferation',        # Proliferating cell population
        'genetic_locus': 'variable_region',     # "variable heavy chain gene"
        'genetic_effect': 'high_variability'    # "high variability" observed
    }

    # Define the properties of each possible answer
    process_properties = {
        'A': {
            'name': 'complement activation',
            'type': 'protein_cascade',
            'genetic_effect': False,
            'reason_incorrect': "Complement activation is a protein-based system and does not cause genetic changes like 'high variability' in B cells."
        },
        'B': {
            'name': 'somatic hypermutation',
            'type': 'genetic_modification',
            'location': 'secondary_lymphoid_organ',
            'timing': 'post_antigen_exposure',
            'genetic_locus': 'variable_region',
            'genetic_effect': 'high_variability',
            'reason_incorrect': "" # This is the correct answer
        },
        'C': {
            'name': 'class switching recombination',
            'type': 'genetic_modification',
            'location': 'secondary_lymphoid_organ',
            'timing': 'post_antigen_exposure',
            'genetic_locus': 'constant_region', # This is the key difference
            'genetic_effect': 'isotype_switch',
            'reason_incorrect': "Class switching recombination affects the *constant* region of the heavy chain gene, not the *variable* region as specified in the question."
        },
        'D': {
            'name': 'VDJ recombination',
            'type': 'genetic_modification',
            'location': 'primary_lymphoid_organ', # Key difference
            'timing': 'pre_antigen_exposure',   # Key difference
            'genetic_locus': 'variable_region',
            'genetic_effect': 'initial_diversity',
            'reason_incorrect': "VDJ recombination occurs in primary lymphoid organs (bone marrow) *before* antigen encounter, not in secondary lymphoid organs *after* antigen exposure as described in the scenario."
        }
    }

    # The correct answer based on our analysis
    correct_answer_letter = 'B'

    if proposed_answer_letter == correct_answer_letter:
        # Check if the proposed answer matches all the facts
        chosen_process = process_properties[proposed_answer_letter]
        
        location_match = chosen_process['location'] == question_facts['location']
        timing_match = chosen_process['timing'] == question_facts['timing']
        locus_match = chosen_process['genetic_locus'] == question_facts['genetic_locus']
        effect_match = chosen_process['genetic_effect'] == question_facts['genetic_effect']

        if location_match and timing_match and locus_match and effect_match:
            return "Correct"
        else:
            # This case should not be reached if our logic is sound, but it's good practice
            return f"The proposed answer {proposed_answer_letter} is identified as correct, but there's a logical inconsistency in the checker. Please review."
    else:
        # The proposed answer is incorrect, provide the specific reason why
        if proposed_answer_letter in process_properties:
            return f"Incorrect. The proposed answer is {proposed_answer_letter} ({process_properties[proposed_answer_letter]['name']}). The reason it is wrong is: {process_properties[proposed_answer_letter]['reason_incorrect']}"
        else:
            return f"Incorrect. The proposed answer letter '{proposed_answer_letter}' is not a valid option."

# --- Provided data from the user ---
question_text = """
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) complement activation
B) somatic hypermutation
C) class switching recombination
D) VDJ recombination
"""

candidate_answers_texts = [
    # ... (all 17 candidate answers would be here) ...
]

final_answer_text = """
Here is a step-by-step analysis of the question to determine the final answer.

**1. Deconstruct the Experimental Scenario**
*   **Antigen and Delivery:** Rotavirus capsid proteins, an antigen, are delivered orally. This route targets the mucosal immune system, specifically the Gut-Associated Lymphoid Tissue (GALT).
*   **Location of Response:** The immune response is observed in the Peyer's patches. These are secondary lymphoid organs embedded in the intestinal wall, which are primary sites for initiating immune responses to antigens from the gut.
*   **Cellular Response:** A "proliferating cell population" is isolated. This indicates clonal expansion of lymphocytes that have recognized the antigen. The mention of a "heavy chain gene" confirms these are B cells, and this proliferation occurs in specialized structures called germinal centers.
*   **Key Genetic Finding:** Sequencing reveals "high variability" specifically in the "variable heavy chain gene." The variable region of an antibody is the part that binds to the antigen.

**2. Evaluate the Options Based on the Key Finding**
The central task is to identify the process that generates high variability in the *variable region* of immunoglobulin genes in *proliferating B cells* within a *secondary lymphoid organ* *after* antigen exposure.

*   **A) complement activation:** This is an effector mechanism of the immune system involving a cascade of plasma proteins. It does not involve any genetic modification of B cells and cannot account for the observed genetic variability. Therefore, this option is incorrect.

*   **B) somatic hypermutation (SHM):** This process occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs (like Peyer's patches). SHM's specific function is to introduce a high rate of point mutations into the genes that code for the *variable regions* of both the heavy and light chains. This directly causes the "high variability" observed in the experiment and allows for the selection of B cells with higher affinity for the antigen (affinity maturation). This perfectly matches all aspects of the scenario.

*   **C) class switching recombination (CSR):** This process also occurs in activated B cells within germinal centers, often concurrently with SHM. However, CSR modifies the *constant* region of the heavy chain gene, not the *variable* region. This changes the antibody's class or isotype (e.g., from IgM to IgA), which alters its function but does not change the antigen-binding site or contribute to the variability observed there. Therefore, this option is incorrect.

*   **D) VDJ recombination:** This is the genetic recombination process that creates the initial, vast diversity of B cell receptors. Crucially, this occurs in the bone marrow (a primary lymphoid organ) during the development of B cells, *before* they ever encounter an antigen. The scenario describes a process happening in mature B cells in a secondary lymphoid organ *after* antigen encounter. Therefore, this option is incorrect.

**3. Conclusion**
The experimental scenario describes a classic germinal center reaction. The combination of a proliferating B cell population in a Peyer's patch after antigen stimulation, coupled with the finding of high genetic variability specifically in the variable heavy chain gene, is the definitive hallmark of somatic hypermutation.

<<<B>>>
"""

# Run the check
result = check_answer(question_text, candidate_answers_texts, final_answer_text)
print(result)