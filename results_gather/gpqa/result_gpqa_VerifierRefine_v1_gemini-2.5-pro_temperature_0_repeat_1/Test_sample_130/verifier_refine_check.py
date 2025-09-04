import re

def check_answer():
    """
    This function checks the correctness of the LLM's answer by codifying the chemical principles
    involved in the problem.

    The logic flow is as follows:
    1.  Identify Reactants: Verify the structures of the dienophile (maleic anhydride) and the diene from the problem description.
    2.  Analyze Reaction & Stereochemistry: Confirm the reaction is a Diels-Alder and determine the major product (*exo* or *endo*) based on steric hindrance, which is a known exception to the endo rule in this case.
    3.  Predict Proton Signals: Assign the described NMR signals to the protons in the product structure.
    4.  Analyze NOESY Data: Based on the 3D structure of the major product, predict which protons will be close in space and thus show a NOESY cross-peak.
    5.  Compare with Options: Match the predicted NOESY correlation to the given multiple-choice options.
    6.  Final Verdict: Check if the LLM's answer matches the derived correct option.
    """

    # --- Step 1: Define the problem constraints and the LLM's answer ---
    question = {
        "dienophile_1H_NMR": {"peaks": 1, "shift_ppm": 7},
        "dienophile_13C_NMR": {"peaks": 2, "shifts_ppm": [137, 165]},
        "diene": "1,2,3,4-tetramethyl-1,3-cyclopentadiene",
        "reaction": "Diels-Alder",
        "noesy_condition": "cross-peak present in major product, absent in minor",
        "options": {
            "A": "A 1H doublet at ~1.5 ppm and a 2H singlet at ~3.5 ppm",
            "B": "A 6H singlet at ~1 ppm and a 6H singlet at ~1.7 ppm",
            "C": "A 6H singlet at ~1 ppm and a 1H doublet at ~1.5 ppm",
            "D": "A 6H singlet at ~1.7 ppm and a 2H singlet at ~3.5 ppm"
        }
    }
    llm_final_answer = "A"

    # --- Step 2: Verify Reactant Identification ---
    # Dienophile: A single 1H peak at 7ppm (vinylic) and two 13C peaks (alkene C=C and anhydride C=O)
    # perfectly describes the symmetric molecule, maleic anhydride.
    dienophile = "maleic anhydride"
    # Diene is given.
    diene = question["diene"]

    # --- Step 3: Analyze Reaction and Stereochemistry ---
    # The reaction is a Diels-Alder. It can form 'endo' and 'exo' products.
    # The endo rule (kinetic preference) is often overridden by severe steric hindrance.
    # The diene has methyl groups at C1 and C4. In the 'endo' transition state, these methyls
    # would clash with the anhydride ring.
    # Therefore, the 'exo' product is sterically favored and is the MAJOR product.
    # The 'endo' product is the MINOR product.
    major_product_conformation = "exo"
    minor_product_conformation = "endo"

    # --- Step 4: Predict Proton Signals and Assign them ---
    # Based on typical chemical shifts for this bicyclic system.
    proton_signals = {
        "bridgehead_methyls": {"description": "6H singlet", "shift_ppm": 1.0, "protons": "Me(C1/C4)"},
        "vinylic_methyls": {"description": "6H singlet", "shift_ppm": 1.7, "protons": "Me(C2/C3)"},
        "bridge_proton": {"description": "1H doublet", "shift_ppm": 1.5, "protons": "H(C7-anti)"},
        "anhydride_protons": {"description": "2H singlet", "shift_ppm": 3.5, "protons": "H(C5/C6)"}
    }

    # --- Step 5: Analyze NOESY Data for Each Isomer ---
    # NOESY shows through-space proximity (< 5 Å).
    # In the MAJOR ('exo') product, the anhydride ring is away from the vinylic methyls.
    # The anhydride protons H(C5/C6) are on the same face as the 'anti' bridge proton H(C7-anti).
    major_product_noe = {"anhydride_protons", "bridge_proton"}

    # In the MINOR ('endo') product, the anhydride ring is tucked under the vinylic methyls.
    # The anhydride protons H(C5/C6) are close to the vinylic methyls Me(C2/C3).
    minor_product_noe = {"anhydride_protons", "vinylic_methyls"}

    # The question asks for the NOE present in the MAJOR product.
    # This is the correlation between "anhydride_protons" and "bridge_proton".
    predicted_signals_in_noe = [
        proton_signals["anhydride_protons"],
        proton_signals["bridge_proton"]
    ]

    # --- Step 6: Compare with Options ---
    derived_correct_option = None
    for option_key, option_text in question["options"].items():
        # Check if the option text contains the descriptions of the two predicted signals
        signal1_desc = predicted_signals_in_noe[0]["description"]
        signal1_shift = str(predicted_signals_in_noe[0]["shift_ppm"])
        signal2_desc = predicted_signals_in_noe[1]["description"]
        signal2_shift = str(predicted_signals_in_noe[1]["shift_ppm"])

        if re.search(signal1_desc, option_text) and re.search(signal1_shift, option_text) and \
           re.search(signal2_desc, option_text) and re.search(signal2_shift, option_text):
            derived_correct_option = option_key
            break

    # --- Step 7: Final Verdict ---
    if derived_correct_option is None:
        return "Error in checking logic: Could not match the derived correct NOE to any of the options."

    if llm_final_answer == derived_correct_option:
        return "Correct"
    else:
        # Explain why the LLM's answer is wrong and what the correct answer is.
        reason = f"The provided answer '{llm_final_answer}' is incorrect. The correct answer is '{derived_correct_option}'.\n"
        reason += f"Reasoning: The major product of this reaction is the 'exo' adduct due to steric hindrance from the methyl groups on the diene.\n"
        reason += f"In the 'exo' structure, the anhydride protons (H at C5/C6, a 2H singlet at ~3.5 ppm) are spatially close to one of the methylene bridge protons (H at C7, a 1H doublet at ~1.5 ppm).\n"
        reason += f"This proximity results in a NOESY cross-peak, which corresponds to option {derived_correct_option}.\n"
        
        # Check if the LLM chose the NOE for the minor product
        if llm_final_answer == 'D':
            reason += "The answer 'D' describes the NOE correlation expected for the MINOR ('endo') product, where the anhydride protons are close to the vinylic methyls (6H singlet at ~1.7 ppm)."
        
        return reason

# Execute the check and print the result
result = check_answer()
print(result)