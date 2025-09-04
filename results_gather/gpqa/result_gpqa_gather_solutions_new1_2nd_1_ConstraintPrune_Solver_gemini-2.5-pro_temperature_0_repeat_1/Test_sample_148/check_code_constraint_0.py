def check_answer_correctness():
    """
    Checks the correctness of the answer to the analytical chemistry question.

    The function models the properties of each answer choice and compares them
    against the constraints derived from the experimental data (NMR, LC, MS).
    """
    # The final answer provided by the LLM to be checked.
    llm_answer = "D"

    # Define the properties of each possible explanation based on chemical principles.
    options = {
        "A": {
            "name": "Mixture of enantiomers",
            "is_isomer": True,  # Same mass as the target molecule.
            "is_separable_by_achiral_lc": False, # Identical physical properties in an achiral environment.
            "is_distinguishable_by_achiral_nmr": False # Identical chemical environment in an achiral solvent.
        },
        "B": {
            "name": "Contaminated with a precursor",
            "is_isomer": False, # Different (lower) mass.
            "is_separable_by_achiral_lc": True, # Different molecule, so different properties.
            "is_distinguishable_by_achiral_nmr": True # Different molecule, so different spectrum.
        },
        "C": {
            "name": "'Double coupling' product",
            "is_isomer": False, # Different (higher) mass.
            "is_separable_by_achiral_lc": True, # Different molecule, so different properties.
            "is_distinguishable_by_achiral_nmr": True # Different molecule, so different spectrum.
        },
        "D": {
            "name": "Mixture of diastereoisomers",
            "is_isomer": True, # Same mass as the target molecule.
            "is_separable_by_achiral_lc": True, # Different physical properties.
            "is_distinguishable_by_achiral_nmr": True # Different chemical environments.
        }
    }

    # Define the constraints derived from the experimental observations in the question.
    experimental_constraints = {
        "is_isomer": True,  # From MS: "Both peaks have the same mass spectrum...consistent with the expected molecule."
        "is_separable_by_achiral_lc": True, # From LC: "two clearly defined peaks."
        "is_distinguishable_by_achiral_nmr": True # From NMR: "two peaks that both correspond to the same alpha-proton."
    }

    # Find the option that satisfies all experimental constraints.
    correct_option = None
    for option_key, properties in options.items():
        if (properties["is_isomer"] == experimental_constraints["is_isomer"] and
            properties["is_separable_by_achiral_lc"] == experimental_constraints["is_separable_by_achiral_lc"] and
            properties["is_distinguishable_by_achiral_nmr"] == experimental_constraints["is_distinguishable_by_achiral_nmr"]):
            # This option's properties match all observations.
            correct_option = option_key
            break # Assuming only one option can be correct.

    # Check if the LLM's answer matches the logically derived correct option.
    if llm_answer == correct_option:
        return "Correct"
    else:
        # Construct a detailed reason for the mismatch.
        reason = f"The provided answer '{llm_answer}' is incorrect.\n"
        reason += f"The correct answer should be '{correct_option}'.\n\n"
        reason += "Analysis based on experimental constraints:\n"
        reason += f"1. MS data (same mass) implies the species are isomers. This rules out options B (precursor) and C (double coupling product).\n"
        reason += f"2. LC data (two peaks) and NMR data (two signals) imply the species are separable and distinguishable in an achiral environment. This rules out option A (enantiomers).\n\n"
        reason += f"Evaluating the options:\n"
        for key, props in options.items():
            if key == correct_option:
                reason += f"- Option {key} ({props['name']}): Matches all constraints.\n"
            else:
                mismatches = []
                if props["is_isomer"] != experimental_constraints["is_isomer"]:
                    mismatches.append("fails the isomer (mass) constraint")
                if props["is_separable_by_achiral_lc"] != experimental_constraints["is_separable_by_achiral_lc"]:
                    mismatches.append("fails the LC separability constraint")
                if props["is_distinguishable_by_achiral_nmr"] != experimental_constraints["is_distinguishable_by_achiral_nmr"]:
                    mismatches.append("fails the NMR distinguishability constraint")
                reason += f"- Option {key} ({props['name']}): Fails because it {', '.join(mismatches)}.\n"

        reason += f"\nTherefore, the only option that satisfies all conditions is '{correct_option}'. The provided answer was '{llm_answer}'."
        return reason

# Execute the check and print the result.
result = check_answer_correctness()
print(result)