def check_spectroscopy_answer():
    """
    Checks the correctness of the identified compound based on spectral data.
    """

    # --- The final answer provided by the LLM to be checked ---
    llm_final_answer = "D"

    # --- Define properties of each candidate compound ---
    # 'nmr_ch2_env' is critical: 'O-CH2' (ester/ether) vs 'CO-CH2' (propionate)
    # 'functional_group_primary' is the main nitrogen-containing group.
    # 'functional_group_carbonyl' is the main carbonyl-containing group.
    compounds = {
        "A": {
            "name": "4-aminophenyl propionate",
            "functional_group_primary": "primary_amine",
            "functional_group_carbonyl": "ester",
            "substitution": "para",
            "nmr_ch2_env": "CO-CH2"  # CH2 is next to a carbonyl, shift ~2.5 ppm
        },
        "B": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "functional_group_primary": "secondary_amide",
            "functional_group_carbonyl": "amide",
            "substitution": "para",
            "nmr_ch2_env": "O-CH2"
        },
        "C": {
            "name": "3-ethoxybenzamide",
            "functional_group_primary": "primary_amide",
            "functional_group_carbonyl": "amide",
            "substitution": "meta",
            "nmr_ch2_env": "O-CH2"
        },
        "D": {
            "name": "ethyl 4-aminobenzoate",
            "functional_group_primary": "primary_amine",
            "functional_group_carbonyl": "ester",
            "substitution": "para",
            "nmr_ch2_env": "O-CH2"  # CH2 is next to an ester oxygen, shift ~4.5 ppm
        }
    }

    # --- Logic to find the correct compound based on spectral data ---
    correct_compound_key = None
    failure_reasons = {}

    for key, props in compounds.items():
        # 1. IR Check: Must be a primary amine (or primary amide for N-H, but check C=O next)
        # The data shows two N-H stretches, ruling out secondary amides.
        if props["functional_group_primary"] not in ["primary_amine", "primary_amide"]:
            failure_reasons[key] = f"Fails IR check: Data indicates a primary amine (2 N-H stretches), but this is a {props['functional_group_primary']}."
            continue

        # 2. IR Check: Must be a conjugated ester (C=O at 1720 cm-1)
        # This rules out amides, which have a C=O stretch < 1700 cm-1.
        if props["functional_group_carbonyl"] != "ester":
            failure_reasons[key] = f"Fails IR check: Data indicates an ester (C=O at 1720 cm-1), but this is an {props['functional_group_carbonyl']}."
            continue

        # 3. NMR Check: Must be para-substituted (two doublets in aromatic region)
        if props["substitution"] != "para":
            failure_reasons[key] = f"Fails NMR check: Data indicates a para-substituted ring, but this is {props['substitution']}-substituted."
            continue

        # 4. NMR Check: Quartet at 4.5 ppm means the ethyl's CH2 is attached to an oxygen.
        if props["nmr_ch2_env"] != "O-CH2":
            failure_reasons[key] = f"Fails NMR check: The quartet at 4.5 ppm indicates an -O-CH2- group. This compound has a {props['nmr_ch2_env']} group, which would have a much lower chemical shift."
            continue

        # If all checks pass, this is the correct compound
        correct_compound_key = key
        break

    # --- Final Verdict ---
    if correct_compound_key is None:
        return "Error: The checker could not identify a matching compound. The spectral data may be inconsistent with the options."

    if correct_compound_key == llm_final_answer:
        return "Correct"
    else:
        reason = failure_reasons.get(llm_final_answer, "it does not satisfy all spectral constraints.")
        return (f"Incorrect. The provided answer is {llm_final_answer}, but the spectral data "
                f"unambiguously points to option {correct_compound_key} ({compounds[correct_compound_key]['name']}).\n"
                f"Reason: The provided answer {llm_final_answer} ({compounds[llm_final_answer]['name']}) is wrong because {reason}")

# Run the checker and print the result
result = check_spectroscopy_answer()
print(result)