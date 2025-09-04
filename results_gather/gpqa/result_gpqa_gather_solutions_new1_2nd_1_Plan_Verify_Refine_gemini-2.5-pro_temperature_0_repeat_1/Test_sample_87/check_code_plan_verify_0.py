import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of an LLM's answer to a multiple-choice immunology question.

    The function models the key characteristics of each immunological process and compares them
    against the constraints derived from the question's experimental scenario.
    """

    # Define the properties of each immunological process based on the options.
    # The options are:
    # A) class switching recombination
    # B) VDJ recombination
    # C) complement activation
    # D) somatic hypermutation
    processes = {
        "A": {
            "name": "class switching recombination",
            "location": "secondary lymphoid organ",
            "timing": "post-antigen",
            "cellular_process": True,
            "genetic_target": "constant region",
            "effect": "isotype change"
        },
        "B": {
            "name": "VDJ recombination",
            "location": "primary lymphoid organ",
            "timing": "pre-antigen",
            "cellular_process": True,
            "genetic_target": "variable region",
            "effect": "initial diversity"
        },
        "C": {
            "name": "complement activation",
            "location": "blood/tissue",
            "timing": "post-antigen",
            "cellular_process": False, # It's a protein cascade, not a genetic event in a cell
            "genetic_target": "none",
            "effect": "pathogen lysis/opsonization"
        },
        "D": {
            "name": "somatic hypermutation",
            "location": "secondary lymphoid organ",
            "timing": "post-antigen",
            "cellular_process": True,
            "genetic_target": "variable region",
            "effect": "high variability"
        }
    }

    # Define the constraints derived from the question's text.
    # "Peyer patches" -> secondary lymphoid organ
    # "isolate the proliferating cell population" -> post-antigen cellular process
    # "sequence their variable heavy chain gene" -> genetic target is the variable region
    # "observe high variability" -> the key effect
    constraints = {
        "location": {
            "value": "secondary lymphoid organ",
            "reason": "The process occurs in Peyer's patches, which are secondary lymphoid organs."
        },
        "timing": {
            "value": "post-antigen",
            "reason": "The process is a response to antigen stimulation in mature lymphocytes."
        },
        "cellular_process": {
            "value": True,
            "reason": "The process involves a 'proliferating cell population' undergoing genetic changes."
        },
        "genetic_target": {
            "value": "variable region",
            "reason": "The question explicitly states high variability was observed in the 'variable heavy chain gene'."
        },
        "effect": {
            "value": "high variability",
            "reason": "The key finding was 'high variability' in the gene sequence."
        }
    }

    # Extract the final answer letter from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>>, <<<B>>>, <<<C>>>, or <<<D>>> in the provided text."

    answer_key = match.group(1)
    selected_process = processes.get(answer_key)

    if not selected_process:
        return f"The answer key '{answer_key}' is invalid."

    # Check the selected process against all constraints
    for key, constraint_info in constraints.items():
        expected_value = constraint_info["value"]
        actual_value = selected_process.get(key)
        
        if actual_value != expected_value:
            reason = constraint_info["reason"]
            return (f"Incorrect. The chosen answer, '{selected_process['name']}', does not satisfy a key constraint. "
                    f"Constraint not met: {reason} "
                    f"The chosen process has '{key}' as '{actual_value}', but the question requires it to be '{expected_value}'.")

    return "Correct"

# The user-provided question and the LLM's answer to be checked.
question = """
You tag the inner capsid proteins of rotavirus with a GFP and deliver them to mice via oral gavage. By in vivo imaging, you observe that some of the fluorescence was destroyed by the stomach acid, and some were entangled in the intestinal mucus before they reached the secondary lymphoid organs. You identify the fluorescent Peyer patches and isolate the proliferating cell population. You sequence their variable heavy chain gene and observe high variability. Which process is taking place?

A) class switching recombination
B) VDJ recombination
C) complement activation
D) somatic hypermutation
"""

llm_answer = """
Here is a step-by-step analysis of the question, incorporating careful points to distinguish between the possible answers.

**1. Deconstruct the Experimental Scenario**
*   **Antigen and Delivery:** Rotavirus proteins are delivered orally. This targets the mucosal immune system, specifically the Gut-Associated Lymphoid Tissue (GALT).
*   **Location of Response:** The immune response is observed in the Peyer's patches. These are secondary lymphoid organs, which are the sites where adaptive immune responses are initiated and mature.
*   **Cellular Response:** A "proliferating cell population" is isolated. This indicates clonal expansion of lymphocytes that have been activated by the antigen. This proliferation occurs in specialized microenvironments called germinal centers.
*   **Genetic Analysis:** The "variable heavy chain gene" is sequenced. The term "heavy chain" confirms the proliferating cells are B lymphocytes.
*   **Key Finding:** The sequencing reveals "high variability" in the *variable* region of this gene. This is the most critical piece of evidence, indicating that new mutations are being actively introduced into the part of the gene that codes for the antigen-binding site.

**2. Analyze the Options with Careful Points**

*   **A) class switching recombination:** This process occurs in activated B cells in secondary lymphoid organs, which matches the location and timing.
    *   **Careful Point:** However, class switching modifies the **constant region** of the heavy chain gene, not the **variable region**. It changes the antibody's function (e.g., from IgM to IgA) but does not alter the antigen-binding site or create the variability observed in the variable region. Therefore, it is incorrect.

*   **B) VDJ recombination:** This is the process that creates the initial diversity of B cell receptors.
    *   **Careful Point:** This occurs during B cell development in a **primary lymphoid organ** (the bone marrow), **before** the B cell ever encounters an antigen. The experiment describes a process happening in a **secondary lymphoid organ** (Peyer's patch) **after** antigen stimulation. Therefore, it is incorrect.

*   **C) complement activation:** This is a system of proteins in the blood that helps clear pathogens.
    *   **Careful Point:** This is a component of the humoral immune system but does not involve any genetic modification within B cells. It cannot explain the observed genetic variability in the heavy chain gene. Therefore, it is incorrect.

*   **D) somatic hypermutation (SHM):** This is a cellular process that introduces point mutations into immunoglobulin genes.
    *   **Careful Point 1 (Location & Timing):** SHM occurs in activated, proliferating B cells within the germinal centers of **secondary lymphoid organs** (like Peyer's patches) **after** antigen stimulation. This perfectly matches the experimental context.
    *   **Careful Point 2 (Genetic Target):** SHM's specific function is to introduce a high rate of mutations into the **variable regions** of the heavy and light chain genes. This is the direct cause of the "high variability" observed in the experiment and is the basis for affinity maturation. This option perfectly matches all aspects of the scenario.

**Conclusion**

The combination of location (Peyer's patch), timing (post-antigen stimulation), cellular activity (proliferation), and specific genetic outcome (high variability in the variable heavy chain gene) are the definitive hallmarks of somatic hypermutation.

<<<D>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)