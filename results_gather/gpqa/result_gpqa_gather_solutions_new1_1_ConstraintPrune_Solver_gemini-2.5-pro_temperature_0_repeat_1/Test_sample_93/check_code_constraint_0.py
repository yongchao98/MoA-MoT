def check_synthesis_answer():
    """
    Checks the correctness of the starting material for the synthesis of
    methyl 2-propyl-1,2,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylate.
    """

    # --- Define Target Molecule Properties ---
    # Based on IUPAC name: methyl 2-propyl-1,2,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylate
    target_properties = {
        "skeleton": "fused_bicyclo_4.4.0",  # decalin
        "saturation": "octahydro",  # one double bond
        "double_bond_position": "C3=C4",
        "regiochemistry": "COOMe_at_C1_next_to_bridgehead_and_propyl_at_C2_next_to_double_bond"
    }

    # --- Define Starting Material Options and Their Predicted Products ---
    options_analysis = {
        "A": {
            "name": "methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate",
            "reaction_type": "IMDA",
            "predicted_product": {
                "skeleton": "fused_bicyclo_4.4.0",
                "saturation": "octahydro",
                "double_bond_position": "C3=C4",
                # Regiochemistry is swapped compared to target
                "regiochemistry": "COOMe_at_C2_next_to_double_bond_and_propyl_at_C1_next_to_bridgehead",
                "reason_for_failure": "produces the wrong regiochemistry (the substituents are swapped relative to the double bond and bridgehead)."
            }
        },
        "B": {
            "name": "methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate",
            "reaction_type": "IMDA",
            "predicted_product": {
                "skeleton": "fused_bicyclo_4.4.0",
                "saturation": "octahydro",
                "double_bond_position": "C3=C4", # Maps from C9=C10 in precursor
                "regiochemistry": "COOMe_at_C1_next_to_bridgehead_and_propyl_at_C2_next_to_double_bond"
            }
        },
        "C": {
            "name": "1-vinylcyclohex-1-ene and methyl hex-2-ynoate",
            "reaction_type": "intermolecular",
            "predicted_product": {
                "skeleton": "fused_bicyclo_4.4.0",
                "saturation": "dihydronaphthalene", # two double bonds
                "reason_for_failure": "uses an alkyne dienophile, producing a product with two double bonds instead of one."
            }
        },
        "D": {
            "name": "Cyclohexene and methyl 2,3-dimethylenehexanoate",
            "reaction_type": "intermolecular",
            "predicted_product": {
                "skeleton": "spirocycle",
                "reason_for_failure": "produces a spirocyclic product, not the required fused bicyclic skeleton."
            }
        }
    }

    # The answer provided by the LLM to be checked
    llm_answer = "B"

    # --- Verification Logic ---
    correct_option = None
    for option, analysis in options_analysis.items():
        product = analysis["predicted_product"]
        is_match = (
            product.get("skeleton") == target_properties["skeleton"] and
            product.get("saturation") == target_properties["saturation"] and
            product.get("double_bond_position") == target_properties["double_bond_position"] and
            product.get("regiochemistry") == target_properties["regiochemistry"]
        )
        if is_match:
            correct_option = option
            break

    if llm_answer == correct_option:
        return "Correct"
    else:
        if llm_answer in options_analysis:
            reason = options_analysis[llm_answer]["predicted_product"].get("reason_for_failure", "does not match the target structure.")
            return f"Incorrect. The provided answer '{llm_answer}' is wrong because it {reason}"
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

# Execute the check and print the result
result = check_synthesis_answer()
print(result)