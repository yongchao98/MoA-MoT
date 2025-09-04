def check_chemistry_answer():
    """
    This function checks the correctness of the answer to the peptide analysis question.
    It simulates the logical deduction process based on the provided analytical data.
    """

    # 1. Define the experimental evidence from the problem description.
    # These are the facts we must satisfy.
    evidence = {
        "species_have_same_mass": True,  # From MS: "Both peaks have the same mass spectrum"
        "species_are_separable_by_lc": True, # From LC: "two clearly defined peaks"
        "species_are_distinguishable_by_nmr": True # From NMR: "two peaks that both correspond to the same alpha-proton"
    }

    # 2. Define the properties of each possible explanation (the options).
    # This represents our knowledge of chemistry. We assume standard (achiral) analytical conditions
    # as none other were specified.
    option_properties = {
        'A': { # Contaminated with a precursor
            "description": "Contamination with a precursor",
            "has_same_mass": False, # A precursor is a different molecule with a different mass.
            "separable_by_lc": True,
            "distinguishable_by_nmr": True
        },
        'B': { # Mixture of enantiomers
            "description": "Mixture of enantiomers",
            "has_same_mass": True, # Enantiomers are isomers, so they have the same mass.
            "separable_by_lc": False, # Enantiomers are not separable by standard (achiral) LC.
            "distinguishable_by_nmr": False # Enantiomers are not distinguishable by standard (achiral) NMR.
        },
        'C': { # 'Double coupling'
            "description": "'Double coupling' has occurred",
            "has_same_mass": False, # A double-coupled product has a higher mass.
            "separable_by_lc": True,
            "distinguishable_by_nmr": True
        },
        'D': { # Mixture of diastereoisomers
            "description": "Mixture of diastereoisomers",
            "has_same_mass": True, # Diastereoisomers are isomers, so they have the same mass.
            "separable_by_lc": True, # Diastereoisomers have different physical properties and are separable.
            "distinguishable_by_nmr": True # Diastereoisomers are chemically distinct and are distinguishable.
        }
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = 'D'

    # 3. Systematically check each option against the evidence.
    valid_options = []
    for option_key, properties in option_properties.items():
        # Check 1: Does the mass property match the MS evidence?
        mass_check = (properties["has_same_mass"] == evidence["species_have_same_mass"])
        
        # Check 2: Does the LC separability match the LC evidence?
        lc_check = (properties["separable_by_lc"] == evidence["species_are_separable_by_lc"])
        
        # Check 3: Does the NMR distinguishability match the NMR evidence?
        nmr_check = (properties["distinguishable_by_nmr"] == evidence["species_are_distinguishable_by_nmr"])

        if mass_check and lc_check and nmr_check:
            valid_options.append(option_key)

    # 4. Determine the final result.
    if not valid_options:
        return "Logic Error: No option satisfies all the experimental evidence."
    
    if len(valid_options) > 1:
        return f"Logic Error: The evidence is consistent with multiple options: {valid_options}."

    derived_correct_answer = valid_options[0]

    if llm_answer == derived_correct_answer:
        return "Correct"
    else:
        # Explain why the LLM's answer is wrong by finding the first failed check.
        llm_properties = option_properties[llm_answer]
        if not (llm_properties["has_same_mass"] == evidence["species_have_same_mass"]):
            reason = f"The chosen answer '{llm_answer}' implies the species have different masses, but the MS data shows they have the same mass."
        elif not (llm_properties["separable_by_lc"] == evidence["species_are_separable_by_lc"]):
            reason = f"The chosen answer '{llm_answer}' implies the species are not separable by standard LC, but the LC data shows two distinct peaks."
        elif not (llm_properties["distinguishable_by_nmr"] == evidence["species_are_distinguishable_by_nmr"]):
            reason = f"The chosen answer '{llm_answer}' implies the species are not distinguishable by standard NMR, but the NMR data shows two distinct signals."
        else:
            reason = "An unknown error occurred in validation."
            
        return f"Incorrect. The provided answer is {llm_answer}, but the correct answer based on the evidence is {derived_correct_answer}. Reason: {reason}"

# Run the check and print the result.
result = check_chemistry_answer()
print(result)