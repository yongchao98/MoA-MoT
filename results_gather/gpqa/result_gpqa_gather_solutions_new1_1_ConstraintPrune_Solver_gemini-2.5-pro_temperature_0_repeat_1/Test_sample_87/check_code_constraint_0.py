import re

def check_immunology_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer by verifying the properties of the selected
    immunological process against the constraints given in the question.

    Args:
        llm_answer_text: The full text of the LLM's response, which should contain an answer
                         in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or the reason for the incorrectness.
    """

    # Define the key constraints derived from the question's text
    constraints = {
        'location': 'secondary lymphoid organ',  # Peyer's patches
        'timing': 'post-antigen encounter',     # Response to antigen delivery
        'genetic_target': 'variable region'     # "variable heavy chain gene"
    }

    # Define the properties of each possible answer choice
    processes = {
        'A': {
            'name': 'VDJ recombination',
            'location': 'primary lymphoid organ',
            'timing': 'pre-antigen encounter',
            'genetic_target': 'variable region',
            'reason_if_wrong': "VDJ recombination occurs in primary lymphoid organs (e.g., bone marrow) before antigen encounter, not in a secondary lymphoid organ as a response to an antigen."
        },
        'B': {
            'name': 'class switching recombination',
            'location': 'secondary lymphoid organ',
            'timing': 'post-antigen encounter',
            'genetic_target': 'constant region',
            'reason_if_wrong': "Class switching recombination affects the constant region of the heavy chain, but the question specifies high variability in the variable region."
        },
        'C': {
            'name': 'somatic hypermutation',
            'location': 'secondary lymphoid organ',
            'timing': 'post-antigen encounter',
            'genetic_target': 'variable region',
            'reason_if_wrong': ""  # This is the correct answer
        },
        'D': {
            'name': 'complement activation',
            'location': 'extracellular',
            'timing': 'post-antigen encounter',
            'genetic_target': 'none',
            'reason_if_wrong': "Complement activation is a protein-based cascade and does not involve genetic mutation in B cells."
        }
    }
    
    # The correct answer key based on immunological facts
    correct_key = 'C'

    # Extract the answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find a valid answer in the format <<<X>>> in the provided text."

    llm_choice_key = match.group(1)
    
    if llm_choice_key != correct_key:
        incorrect_process = processes[llm_choice_key]
        return f"Incorrect. The provided answer is {llm_choice_key} ({incorrect_process['name']}). {incorrect_process['reason_if_wrong']}"

    # Verify that the correct answer indeed satisfies all constraints
    correct_process = processes[correct_key]
    for key, value in constraints.items():
        if correct_process[key] != value:
            # This is an internal sanity check for the script's logic
            return f"Error in checker logic: The designated correct answer '{correct_key}' does not satisfy the constraint '{key}'. Expected '{value}', but got '{correct_process[key]}'."

    return "Correct"

# The final answer provided by the LLM to be checked
final_llm_answer = """
Here is a step-by-step analysis of the question to determine the correct process.

**Step 1: Deconstruct the experimental scenario and key observations.**
*   **Antigen and Delivery:** The experiment uses rotavirus capsid proteins (an antigen) delivered orally. This targets the mucosal immune system of the gastrointestinal tract.
*   **Location of Response:** The immune response is observed in the Peyer's patches. These are secondary lymphoid organs in the small intestine, which are primary sites for initiating adaptive immune responses to gut antigens.
*   **Cellular Response:** A "proliferating cell population" is isolated. This indicates that lymphocytes (in this case, B cells, as implied by the "heavy chain gene") have recognized the antigen and are undergoing clonal expansion.
*   **Genetic Observation:** The "variable heavy chain gene" of these proliferating cells is sequenced and shows "high variability." This is the most critical piece of information. It means that within the population of expanding cells, many new mutations have been introduced into the specific part of the gene that codes for the antigen-binding site.

**Step 2: Evaluate the candidate processes against the observations.**
The question asks to identify the process responsible for generating high variability in the variable region of the heavy chain gene *after* antigen stimulation in a secondary lymphoid organ.

*   **A) VDJ recombination:** This process creates the initial diversity of B cell receptors.
    *   **Careful Point:** This occurs during B cell development in the bone marrow, a *primary* lymphoid organ, *before* the B cell ever encounters an antigen. The experiment is observing events in a *secondary* lymphoid organ *after* antigen encounter. Therefore, VDJ recombination is incorrect.

*   **B) Class switching recombination (CSR):** This process changes the type of antibody produced (e.g., from IgM to IgA).
    *   **Careful Point:** CSR occurs in activated B cells in secondary lymphoid organs, which matches the location and timing. However, it modifies the *constant* region of the heavy chain gene, not the *variable* region. The observation was high variability in the *variable* region. Therefore, CSR is not the direct cause of the observed phenomenon.

*   **C) Somatic hypermutation (SHM):** This is a process that introduces point mutations into immunoglobulin genes.
    *   **Careful Point:** SHM occurs in activated, proliferating B cells within germinal centers of secondary lymphoid organs (like Peyer's patches). Its specific function is to introduce mutations into the *variable* regions of the heavy and light chain genes. This process is the direct cause of the "high variability" observed and is essential for affinity maturation, where antibodies with higher binding strength to the antigen are selected. This matches all aspects of the experimental scenario perfectly.

*   **D) Complement activation:** This is a system of proteins in the blood that helps clear pathogens.
    *   **Careful Point:** This is a component of the humoral immune system but does not involve any genetic modification of B cells. It cannot explain the observed genetic variability. Therefore, complement activation is incorrect.

**Step 3: Synthesize and conclude.**
The combination of location (Peyer's patch), timing (post-antigen stimulation), cellular activity (proliferation), and genetic outcome (high variability in the variable heavy chain gene) are the classic hallmarks of somatic hypermutation.

<<<C>>>
"""

# Run the check
result = check_immunology_answer(final_llm_answer)
print(result)