import re

def check_answer_correctness(question_details, final_answer_text):
    """
    Checks the correctness of the final answer for the given immunology question.

    Args:
        question_details (dict): A dictionary containing the key facts from the question.
        final_answer_text (str): The full text of the final answer provided, including the rationale and the <<<answer>>> tag.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error if it's incorrect.
    """

    # Define the correct logic based on established immunology principles.
    # The question describes a process with the following characteristics:
    # 1. Location: Secondary lymphoid organ (Peyer's patches)
    # 2. Timing: Post-antigen exposure (in proliferating cells)
    # 3. Genetic Target: Variable region of the heavy chain gene
    # 4. Outcome: High variability (point mutations)
    # This combination of features uniquely identifies Somatic Hypermutation (SHM).

    process_definitions = {
        "somatic hypermutation": {
            "option": "A",
            "matches_question": True,
            "reason": "This process occurs in proliferating B cells in secondary lymphoid organs after antigen exposure, introducing high variability into the variable region of immunoglobulin genes."
        },
        "class switching recombination": {
            "option": "B",
            "matches_question": False,
            "reason": "This process affects the *constant* region of the heavy chain, not the *variable* region where high variability was observed."
        },
        "VDJ recombination": {
            "option": "C",
            "matches_question": False,
            "reason": "This process occurs during B cell development in the bone marrow *before* antigen encounter, not in proliferating cells after antigen exposure."
        },
        "complement activation": {
            "option": "D",
            "matches_question": False,
            "reason": "This is a protein-based system and does not involve genetic mutation in B cells."
        }
    }
    
    # The correct option based on the question's details is 'A'.
    correct_option = "A"

    # Extract the chosen answer from the final answer text
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a valid answer in the format <<<A>>>, <<<B>>>, etc., in the provided text."

    chosen_option = match.group(1)

    # Check if the chosen option is the correct one
    if chosen_option == correct_option:
        # Additionally, verify that the rationale provided in the text supports the correct answer.
        # A simple check is to see if the text correctly dismisses other options or correctly describes the chosen one.
        rationale_check_shm = "variable region" in final_answer_text and "high variability" in final_answer_text
        rationale_check_csr = "constant region" in final_answer_text
        rationale_check_vdj = "before antigen encounter" in final_answer_text or "bone marrow" in final_answer_text
        
        if rationale_check_shm and rationale_check_csr and rationale_check_vdj:
             return "Correct"
        else:
            return f"The chosen option {chosen_option} is correct, but the provided rationale seems incomplete or flawed."

    else:
        # Find the reason why the chosen option is incorrect.
        for process_name, details in process_definitions.items():
            if details["option"] == chosen_option:
                return f"Incorrect. The final answer chose option {chosen_option} ({process_name}), but the correct answer is {correct_option} (somatic hypermutation). The reason the chosen answer is wrong is: {details['reason']}"
        
        return f"Incorrect. The final answer chose option {chosen_option}, but the correct answer is {correct_option} (somatic hypermutation)."


# --- Input Data ---

# Key details from the question prompt
question_info = {
    "task": "Decide the final answer based on the following answers from other agents.",
    "question": "You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?",
    "options": {
        "A": "somatic hypermutation",
        "B": "class switching recombination",
        "C": "VDJ recombination",
        "D": "complement activation"
    }
}

# The final answer provided by the user to be checked
final_answer_from_prompt = """
Here is a step-by-step analysis of the question to determine the correct process.

**1. Deconstruct the Experimental Scenario:**
*   **Antigen and Delivery:** Rotavirus capsid proteins are delivered orally. This route targets the mucosal immune system, specifically the Gut-Associated Lymphoid Tissue (GALT).
*   **Location of Immune Response:** The response is observed in the Peyer's patches. These are secondary lymphoid organs embedded in the intestinal wall, which are primary sites for initiating immune responses to antigens from the gut.
*   **Cellular Activity:** A "proliferating cell population" is isolated. This indicates clonal expansion of lymphocytes that have recognized the antigen. In secondary lymphoid organs, this proliferation occurs in specialized structures called germinal centers.
*   **Genetic Analysis:** The "variable heavy chain gene" is sequenced. The term "heavy chain" is specific to B cells, which produce B cell receptors and antibodies. The "variable" region is the part of the antibody that binds to the antigen.
*   **Key Finding:** The sequencing reveals "high variability" in this gene. This means that among the population of proliferating B cells, numerous point mutations have been introduced into the gene segment that determines antigen binding.

**2. Evaluate the Options Based on the Evidence:**
The central task is to identify the process that generates high variability in the *variable region* of immunoglobulin genes in a *proliferating* B cell population within a *secondary lymphoid organ* *after* antigen exposure.

*   **A) Somatic hypermutation (SHM):** This process occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs (like Peyer's patches). It deliberately introduces a high rate of point mutations specifically into the genes encoding the variable regions of antibodies. This process is the direct cause of the "high variability" observed and is essential for affinity maturation, where antibodies with higher binding strength are selected. This perfectly matches all aspects of the experimental scenario.

*   **B) Class switching recombination (CSR):** This process also occurs in activated B cells in germinal centers. However, CSR modifies the *constant* region of the heavy chain gene, not the *variable* region. This changes the antibody's function (e.g., switching from IgM to IgA, which is common in mucosal responses) but does not alter its antigen-binding site or explain the variability found there.

*   **C) VDJ recombination:** This is the genetic recombination process that creates the initial diversity of B cell receptors. Crucially, this happens during the development of B cells in the bone marrow (a primary lymphoid organ), *before* they ever encounter an antigen. The process described in the question is happening in mature B cells in a secondary lymphoid organ *in response* to an antigen.

*   **D) Complement activation:** This is a system of proteins in the blood and tissues that is part of the innate immune response. It helps to clear pathogens but is not a cellular process and does not involve any genetic changes in B cells.

**3. Conclusion:**
The experimental scenario describes a classic germinal center reaction. The combination of a proliferating B cell population in a Peyer's patch showing high genetic variability specifically in the variable heavy chain gene is the definitive signature of somatic hypermutation.

<<<A>>>
"""

# Execute the check
result = check_answer_correctness(question_info, final_answer_from_prompt)
print(result)