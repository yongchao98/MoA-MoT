import re

def check_answer(question, candidate_answers_text, final_answer_text):
    """
    Checks the correctness of the final answer to a multiple-choice immunology question.

    The function programmatically models the constraints given in the question and the
    properties of the possible answers. It then determines the correct answer based on
    these models and compares it to the provided final answer.

    Args:
        question (str): The question text.
        candidate_answers_text (str): A string containing all candidate answers.
        final_answer_text (str): The final answer provided by the agent.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Step 1: Define the constraints from the question ---
    # The question describes a specific experimental scenario. We can extract key facts
    # that the correct answer must satisfy.
    scenario_constraints = {
        "location": "secondary_lymphoid_organ",  # Response is in Peyer's patches.
        "timing": "post_antigen_exposure",      # Response is in "proliferating cells" after antigen challenge.
        "genetic_target": "variable_region",    # "high variability" is in the "variable heavy chain gene".
        "is_b_cell_genetic_process": True       # The process involves genetic changes in B cells.
    }

    # --- Step 2: Define the properties of each possible answer ---
    # Each option corresponds to a distinct biological process with known characteristics.
    processes = {
        "A": {
            "name": "somatic hypermutation",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_exposure",
            "genetic_target": "variable_region",
            "is_b_cell_genetic_process": True
        },
        "B": {
            "name": "VDJ recombination",
            "location": "primary_lymphoid_organ",   # Mismatch: Occurs in bone marrow.
            "timing": "pre_antigen_exposure",      # Mismatch: Occurs during B cell development.
            "genetic_target": "variable_region",
            "is_b_cell_genetic_process": True
        },
        "C": {
            "name": "class switching recombination",
            "location": "secondary_lymphoid_organ",
            "timing": "post_antigen_exposure",
            "genetic_target": "constant_region",   # Mismatch: Affects the constant, not variable, region.
            "is_b_cell_genetic_process": True
        },
        "D": {
            "name": "complement activation",
            "is_b_cell_genetic_process": False     # Mismatch: Not a genetic process within B cells.
        }
    }

    # --- Step 3: Determine the correct answer based on the model ---
    correct_option = None
    for option, properties in processes.items():
        match = True
        for constraint, required_value in scenario_constraints.items():
            if properties.get(constraint) != required_value:
                match = False
                break
        if match:
            correct_option = option
            break

    # --- Step 4: Extract the LLM's final answer ---
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find the final answer in the expected format '<<<A>>>'."
    llm_answer = match.group(1)

    # --- Step 5: Compare the LLM's answer with the derived correct answer ---
    if llm_answer == correct_option:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness
        reasons = []
        llm_choice_properties = processes.get(llm_answer)
        if not llm_choice_properties:
             return f"The provided answer '{llm_answer}' is not a valid option."

        for constraint, required_value in scenario_constraints.items():
            actual_value = llm_choice_properties.get(constraint)
            if actual_value != required_value:
                if constraint == "location":
                    reasons.append(f"the process must occur in a secondary lymphoid organ (like Peyer's patches), but {llm_choice_properties['name']} occurs in a {actual_value}.")
                elif constraint == "timing":
                    reasons.append(f"the process must occur post-antigen exposure in proliferating cells, but {llm_choice_properties['name']} occurs {actual_value}.")
                elif constraint == "genetic_target":
                    reasons.append(f"the question specifies variability in the 'variable region', but {llm_choice_properties['name']} affects the '{actual_value}'.")
                elif constraint == "is_b_cell_genetic_process":
                    reasons.append(f"{llm_choice_properties['name']} is not a genetic modification process within B cells.")

        reason_string = " and ".join(reasons)
        return (f"The answer '{llm_answer}' is incorrect. The correct answer is '{correct_option}'. "
                f"The chosen answer is wrong because {reason_string}")


# The user-provided data
question = """
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) somatic hypermutation
B) VDJ recombination
C) class switching recombination
D) complement activation
"""

candidate_answers_text = """
Answer 1: ... <<<C>>>
Answer 2: ... <<<A>>>
Answer 3: ... <<<B>>>
Answer 4: ... <<<B>>>
Answer 5: ... <<<A>>>
Answer 6: ... (no answer)
Answer 7: ... <<<A>>>
Answer 8: ... <<<D>>>
Answer 9: ... <<<C>>>
Answer 10: ... <<<B>>>
Answer 11: ... <<<B>>>
Answer 12: ... <<<C>>>
Answer 13: ... <<<B>>>
Answer 14: ... <<<D>>>
Answer 15: ... <<<C>>>
"""

final_answer_text = """
Here is a step-by-step analysis of the task based on the provided information.

### Step 1: Deconstruct the Experimental Scenario and Establish Constraints

The question describes a specific immunology experiment. We can break it down into a series of facts that act as constraints to prune the possible answers.

*   **Constraint 1: Antigen and Delivery Route.** The antigen is a rotavirus protein delivered orally. This targets the mucosal immune system of the gut.
*   **Constraint 2: Location of the Immune Response.** The response is observed in the "Peyer's patches," which are secondary lymphoid organs in the intestine, known for initiating immune responses to gut antigens. They contain germinal centers.
*   **Constraint 3: Cellular Activity.** A "proliferating cell population" is isolated. This indicates clonal expansion, a process where lymphocytes that recognize an antigen multiply rapidly. This happens in germinal centers.
*   **Constraint 4: Genetic Target.** The gene being analyzed is the "variable heavy chain gene." This is the part of the B-cell receptor/antibody gene that codes for the antigen-binding site.
*   **Constraint 5: Key Observation.** The finding is "high variability" in this specific gene among the proliferating cells. The correct process must explain the introduction of new mutations in this population.

### Step 2: Evaluate Each Option Against the Constraints

We will now test each option against the established constraints. The process must occur in proliferating cells within Peyer's patches after antigen exposure and cause high variability in the *variable* region of the heavy chain gene.

**A) somatic hypermutation (SHM)**
*   **Mechanism:** Introduces a high rate of point mutations into the variable regions of immunoglobulin genes.
*   **Location:** Occurs in activated, proliferating B cells within the germinal centers of secondary lymphoid organs, including Peyer's patches.
*   **Timing:** Triggered *after* a B cell encounters its specific antigen.
*   **Outcome:** Directly causes "high variability" in the "variable heavy chain gene" to improve antigen binding affinity.
*   **Verdict:** **Passes.** This option perfectly matches all constraints.

**B) VDJ recombination**
*   **Mechanism:** A process of gene segment rearrangement that creates the initial diversity of B-cell receptors.
*   **Location:** Occurs in the bone marrow (a primary lymphoid organ).
*   **Timing:** Occurs during B-cell development, *before* the B cell ever encounters an antigen.
*   **Outcome:** While it creates variability, it is not the process happening in an already-activated, proliferating cell population in a secondary lymphoid organ.
*   **Verdict:** **Fails.** Fails on location (Constraint 2) and timing/cellular context (Constraint 3).

**C) class switching recombination (CSR)**
*   **Mechanism:** A process that changes the *constant* region of the heavy chain gene, not the variable region.
*   **Location:** Occurs in activated, proliferating B cells in germinal centers, often alongside SHM.
*   **Timing:** Occurs *after* antigen encounter.
*   **Outcome:** Changes the antibody's function (e.g., from IgM to IgA) but does not alter the antigen-binding site or cause variability in the *variable* region.
*   **Verdict:** **Fails.** Fails on the genetic target (Constraint 4) and the key observation (Constraint 5).

**D) complement activation**
*   **Mechanism:** A cascade of plasma proteins that helps clear pathogens. It is not a process of genetic modification.
*   **Location:** Occurs in blood and tissues where pathogens are present.
*   **Timing:** Can be triggered by antibodies, but it is not a process within B cells that alters their genes.
*   **Outcome:** Does not involve genetic changes in lymphocytes.
*   **Verdict:** **Fails.** Fails on the fundamental mechanism; it is unrelated to genetic variability in B cells (Constraints 3, 4, 5).

### Step 3: Conclusion

The experimental scenario describes a classic germinal center reaction in the gut-associated lymphoid tissue (GALT). The key finding—high variability in the variable heavy chain gene of proliferating B cells—is the defining characteristic of somatic hypermutation. This process is essential for affinity maturation, where the immune system refines its antibody response to be more effective against a specific pathogen.

<<<A>>>
"""

# Run the check
result = check_answer(question, candidate_answers_text, final_answer_text)
print(result)