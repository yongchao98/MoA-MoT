import re

def check_diels_alder_product():
    """
    Checks the correctness of the selected product for the reaction of
    5-fluorocyclopenta-1,3-diene with maleic anhydride.
    """
    # The options provided in the question
    options = {
        "A": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
    }

    # The answer provided by the LLM
    llm_answer = "C"

    # --- Define Chemical Constraints ---

    # Constraint 1: The major product must be the 'endo' adduct.
    # For this system, 'endo' corresponds to the (4S, 7R) relative stereochemistry.
    def is_endo(name):
        return '4S' in name and '7R' in name

    # Constraint 2: The major product must result from 'syn' attack.
    # For this system, 'syn' attack with respect to the fluorine atom results in an '8s' configuration.
    def is_syn(name):
        # We search for '8s' specifically, not just 's', to avoid confusion with other stereocenters.
        match = re.search(r'8s', name)
        return match is not None

    # --- Apply Constraints to Find the Major Product ---
    major_product_candidates = []
    for key, name in options.items():
        # A product is the major product if it satisfies both constraints.
        if is_endo(name) and is_syn(name):
            major_product_candidates.append(key)

    # --- Validate the LLM's Answer ---
    # There should be exactly one major product among the options.
    if len(major_product_candidates) != 1:
        return f"Error in logic: Expected to find exactly one major product, but found {len(major_product_candidates)} ({major_product_candidates}). The constraints may be incorrectly defined or the options may be ambiguous."

    predicted_correct_answer = major_product_candidates[0]

    if predicted_correct_answer == llm_answer:
        return "Correct"
    else:
        # Find the properties of the LLM's incorrect choice to explain the error.
        llm_choice_name = options[llm_answer]
        llm_choice_is_endo = "endo" if is_endo(llm_choice_name) else "exo"
        llm_choice_is_syn = "syn (8s)" if is_syn(llm_choice_name) else "anti (8r)"

        return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer based on chemical principles is {predicted_correct_answer}.\n"
                f"Reasoning: The major product must be both 'endo' and 'syn'.\n"
                f"- The 'endo' rule eliminates options A and B.\n"
                f"- The 'syn' rule (favoring 8s configuration) eliminates option D.\n"
                f"- Only option {predicted_correct_answer} is both 'endo' and 'syn'.\n"
                f"The chosen answer {llm_answer} is {llm_choice_is_endo} and {llm_choice_is_syn}, which does not match the expected major product if it were incorrect.")

# Execute the check
result = check_diels_alder_product()
print(result)