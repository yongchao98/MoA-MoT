def check_answer_correctness():
    """
    Checks the correctness of the final answer by logically applying the constraints
    from the question to each possible option.
    """
    # The final answer provided by the LLM being checked.
    # The provided answer is <<<B>>>, so we extract 'B'.
    llm_answer = "B"

    # --- Define Constraints from the Question ---
    # 1. Mass Spec: Two peaks, SAME mass, consistent with expected molecule.
    #    This means the components are ISOMERS.
    constraint_mass = "isomers"  # Must have the same mass as the expected product.

    # 2. LC & NMR: Two distinct peaks in LC and two distinct signals in NMR.
    #    This means the components are distinguishable and separable by standard achiral methods.
    constraint_separability = "distinguishable_by_achiral_methods"

    # --- Define Properties of Each Option ---
    options_properties = {
        "A": {
            "description": "Contamination with a precursor",
            "mass_property": "different_mass",  # A precursor is not an isomer.
            "separability_property": "distinguishable_by_achiral_methods"
        },
        "B": {
            "description": "Mixture of diastereoisomers",
            "mass_property": "isomers",  # Diastereomers are isomers.
            "separability_property": "distinguishable_by_achiral_methods" # Diastereomers have different properties.
        },
        "C": {
            "description": "'Double coupling' side product",
            "mass_property": "different_mass",  # A side product with a different formula is not an isomer.
            "separability_property": "distinguishable_by_achiral_methods"
        },
        "D": {
            "description": "Mixture of enantiomers",
            "mass_property": "isomers",  # Enantiomers are isomers.
            "separability_property": "indistinguishable_by_achiral_methods" # Enantiomers are not separable by achiral LC/NMR.
        }
    }

    # --- Check the LLM's Answer against the Constraints ---
    chosen_option = options_properties.get(llm_answer)

    if not chosen_option:
        return f"Invalid answer option '{llm_answer}'. The options are A, B, C, D."

    # Check Mass Constraint
    if chosen_option["mass_property"] != constraint_mass:
        return (f"The answer '{llm_answer}' ({chosen_option['description']}) is incorrect. "
                f"It violates the mass spectrometry constraint. The MS data shows both peaks have the "
                f"same mass as the expected molecule, meaning they must be isomers. "
                f"A {chosen_option['description']} would have a different mass.")

    # Check Separability/Distinguishability Constraint
    if chosen_option["separability_property"] != constraint_separability:
        return (f"The answer '{llm_answer}' ({chosen_option['description']}) is incorrect. "
                f"It violates the LC and NMR constraints. The data shows two distinct, separable peaks. "
                f"However, a {chosen_option['description']} would be indistinguishable by standard "
                f"(achiral) LC and NMR techniques.")

    # If all constraints are satisfied
    return "Correct"

# Run the check and print the result
result = check_answer_correctness()
print(result)