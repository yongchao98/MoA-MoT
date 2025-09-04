import re

def check_correctness_of_answer():
    """
    This function checks the correctness of an LLM's answer to a multiple-choice immunology question.
    It analyzes the question's premises and compares them against the known characteristics of the biological processes listed in the options.
    """

    # The question and the provided answer from another LLM
    question = {
        "text": "You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?",
        "options": {
            "A": "VDJ recombination",
            "B": "class switching recombination",
            "C": "complement activation",
            "D": "somatic hypermutation"
        }
    }
    
    llm_answer = "Thank you for the confirmation. The provided information supports the conclusion that the process taking place is somatic hypermutation. I am ready for the next question."

    # --- Step 1: Parse the LLM's answer to identify the chosen option ---
    llm_choice = None
    for key, value in question["options"].items():
        # Use regex to find the full name of the process in the answer text
        if re.search(r'\b' + re.escape(value) + r'\b', llm_answer, re.IGNORECASE):
            llm_choice = key
            break
    
    if llm_choice is None:
        return "Failure to parse: The LLM's answer does not explicitly mention one of the provided options."

    # --- Step 2: Analyze the question to determine the correct answer based on scientific principles ---
    
    # Key constraints derived from the question:
    # 1. Location: "Peyer patches" are secondary lymphoid organs where germinal center reactions occur.
    # 2. Cell Type: "proliferating cell population" indicates cells undergoing clonal expansion, a hallmark of the germinal center reaction.
    # 3. Genetic Locus: "variable heavy chain gene" is the specific part of the immunoglobulin gene that encodes the antigen-binding site.
    # 4. Observation: "high variability" among this gene in the proliferating population points to a process that introduces mutations.

    # Evaluate each option against the constraints:
    # - A) VDJ recombination: This process creates initial diversity in developing B-cells in the *primary* lymphoid organs (bone marrow), *before* they encounter an antigen. It does not cause ongoing variability in a *proliferating* population in a *secondary* lymphoid organ. This is incorrect.
    # - B) Class switching recombination: This process occurs in germinal centers but alters the *constant* region of the heavy chain, not the *variable* region. It changes the antibody's function (e.g., IgM to IgG), not its antigen-binding specificity or variability. This is incorrect.
    # - C) Complement activation: This is a system of proteins in the blood and tissues (part of the humoral innate immune system). It is not a genetic process occurring within B-cells. This is incorrect.
    # - D) Somatic hypermutation (SHM): This is the precise biological process where activated B-cells in the germinal centers of secondary lymphoid organs (like Peyer's patches) undergo rapid proliferation and introduce point mutations into their variable region genes. This leads to "high variability" and is the basis for affinity maturation. This is correct.

    correct_choice = 'D'

    # --- Step 3: Compare the LLM's choice with the determined correct answer ---
    if llm_choice == correct_choice:
        return "Correct"
    else:
        reason = f"Incorrect. The LLM identified the process as '{question['options'][llm_choice]}' (Option {llm_choice}), but the correct answer is '{question['options'][correct_choice]}' (Option {correct_choice}).\n"
        reason += "Reasoning: The question describes the generation of 'high variability' in the 'variable heavy chain gene' of 'proliferating cells' within 'Peyer's patches'. "
        reason += "This set of conditions—a secondary lymphoid organ, proliferating B-cells, and mutation in the variable gene region—is the specific definition of somatic hypermutation. "
        if llm_choice == 'A':
            reason += "VDJ recombination occurs before antigen encounter in the bone marrow, not in proliferating cells in a Peyer's patch."
        elif llm_choice == 'B':
            reason += "Class switching affects the constant region of the gene, not the variable region where variability was observed."
        elif llm_choice == 'C':
            reason += "Complement activation is a protein cascade, not a genetic process in B-cells."
        return reason

# The final output is the result of executing the checking code.
# In this case, the LLM's answer correctly identifies "somatic hypermutation" (Option D),
# which aligns with the scientific analysis.
result = check_correctness_of_answer()
print(result)