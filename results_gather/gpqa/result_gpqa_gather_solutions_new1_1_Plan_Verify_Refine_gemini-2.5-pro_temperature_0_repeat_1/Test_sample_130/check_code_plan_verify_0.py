import re

def check_answer():
    """
    Checks the correctness of the LLM's answer by codifying the chemical reasoning.
    """
    # 1. Define proton groups and their expected NMR signatures based on the problem description.
    protons = {
        "H_anhydride": {"desc": "2H singlet at ~3.5 ppm", "role": "Protons on the anhydride ring"},
        "H_bridge": {"desc": "1H doublet at ~1.5 ppm", "role": "A proton on the C7 bridge"},
        "Me_vinyl": {"desc": "6H singlet at ~1.7 ppm", "role": "Vinylic methyl groups"},
        "Me_bridgehead": {"desc": "6H singlet at ~1.0 ppm", "role": "Bridgehead methyl groups"}
    }

    # 2. Define key spatial proximities (< 5 Angstroms) that differ between isomers.
    # This is based on the correct 3D models of the products.
    # In the endo product, the anhydride protons (exo face) are near the C7 bridge anti-proton (exo face).
    endo_specific_interaction = ("H_anhydride", "H_bridge")
    # In the exo product, the anhydride protons (endo face) are near the vinylic methyls (endo face).
    exo_specific_interaction = ("H_anhydride", "Me_vinyl")

    # 3. Apply the Alder-endo rule: The major product is the 'endo' adduct.
    # The question asks for a cross-peak present in the major product but absent in the minor.
    # This corresponds to the interaction specific to the endo isomer.
    expected_interaction_protons = endo_specific_interaction

    # 4. Map the interacting protons to their NMR descriptions.
    proton1_desc = protons[expected_interaction_protons[0]]["desc"]
    proton2_desc = protons[expected_interaction_protons[1]]["desc"]

    # 5. Define the options provided in the question.
    options = {
        "A": {"desc_1": "1H doublet at ~1.5 ppm", "desc_2": "2H singlet at ~3.5 ppm"},
        "B": {"desc_1": "6H singlet at ~1.0 ppm", "desc_2": "1H doublet at ~1.5 ppm"},
        "C": {"desc_1": "6H singlet at ~1.7 ppm", "desc_2": "2H singlet at ~3.5 ppm"},
        "D": {"desc_1": "6H singlet at ~1.0 ppm", "desc_2": "6H singlet at ~1.7 ppm"}
    }

    # 6. Determine the correct option letter based on the analysis.
    correct_option = None
    for option, descriptions in options.items():
        if (descriptions["desc_1"] == proton1_desc and descriptions["desc_2"] == proton2_desc) or \
           (descriptions["desc_1"] == proton2_desc and descriptions["desc_2"] == proton1_desc):
            correct_option = option
            break
    
    # 7. The final answer provided by the LLM being checked.
    llm_answer = "A"

    # 8. Compare the derived correct option with the LLM's answer.
    if llm_answer == correct_option:
        return "Correct"
    else:
        reason = f"The provided answer is {llm_answer}, but the correct answer is {correct_option}.\n"
        reason += "Reasoning:\n"
        reason += "1. The major product of this Diels-Alder reaction is the *endo* adduct, according to the Alder-endo rule.\n"
        reason += "2. The question asks for a NOESY cross-peak present in the major product but absent in the minor. This requires finding a spatial proximity unique to the *endo* isomer.\n"
        reason += f"3. In the *endo* isomer, the {protons[expected_interaction_protons[0]]['role']} ({proton1_desc}) are uniquely close to {protons[expected_interaction_protons[1]]['role']} ({proton2_desc}).\n"
        reason += f"4. In the minor (*exo*) isomer, these protons are far apart. Instead, the anhydride protons are close to the vinylic methyl groups ({protons['Me_vinyl']['desc']}).\n"
        reason += f"5. Therefore, the distinguishing cross-peak connects the resonances described in option {correct_option}."
        
        # Check for common errors
        if llm_answer == "C":
             reason += "\nThe provided answer 'C' incorrectly identifies the interaction for the minor (*exo*) product as the one for the major product."
        
        return reason

# Execute the check and print the result.
result = check_answer()
print(result)