def predict_electrocyclic_product_ratio():
    """
    Predicts the product ratio of a thermal electrocyclic reaction
    using Frontier Molecular Orbital (FMO) theory.
    """

    # 1. Define system parameters
    pi_electrons = 8
    condition = "thermal"
    # (2Z,8E) stereochemistry corresponds to an "in/out" arrangement of the
    # terminal methyl substituents.
    reactant_stereochem = "in/out"
    product_A_isomer = "cis"
    product_B_isomer = "trans"

    # 2. Determine the allowed mode of rotation using Woodward-Hoffmann rules
    # Check if the system is 4n or 4n+2
    is_4n = (pi_electrons % 4 == 0)

    rotation_mode = ""
    if condition == "thermal":
        if is_4n:
            # For 4n systems, thermal reactions are conrotatory
            rotation_mode = "conrotatory"
        else: # 4n+2 system
            rotation_mode = "disrotatory"
    # (Other conditions like photochemical would go here)

    # 3. Predict the stereochemistry of the major product
    predicted_product_isomer = ""
    if rotation_mode == "conrotatory":
        if reactant_stereochem == "in/out":
            predicted_product_isomer = "trans"
        else: # "in/in" or "out/out"
            predicted_product_isomer = "cis"
    elif rotation_mode == "disrotatory":
        if reactant_stereochem == "in/out":
            predicted_product_isomer = "cis"
        else: # "in/in" or "out/out"
            predicted_product_isomer = "trans"

    # 4. Determine the theoretical product ratio
    # FMO theory predicts the "allowed" product is formed exclusively.
    # The "forbidden" product is not formed.
    ratio_A = 0
    ratio_B = 0
    if predicted_product_isomer == product_A_isomer:
        # Product A is the allowed product
        ratio_A = 1
        ratio_B = 0
        allowed_pathway_product = "A (cis)"
        forbidden_pathway_product = "B (trans)"
    else:
        # Product B is the allowed product
        ratio_A = 0
        ratio_B = 1
        allowed_pathway_product = "B (trans)"
        forbidden_pathway_product = "A (cis)"

    # 5. Print the reasoning and the final result
    print("--- FMO Theory Analysis ---")
    print(f"Reactant System: {pi_electrons} π-electrons")
    print(f"Reaction Condition: {condition.capitalize()}")
    print(f"Reactant Stereochemistry: {reactant_stereochem.capitalize()}")
    print("\n--- Prediction ---")
    print(f"Based on Woodward-Hoffmann rules, a {pi_electrons} π-electron system under thermal conditions undergoes a '{rotation_mode}' ring closure.")
    print(f"A '{rotation_mode}' closure of a reactant with '{reactant_stereochem}' stereochemistry yields a '{predicted_product_isomer}' product.")
    print(f"The 'allowed' pathway leads to the {allowed_pathway_product} isomer.")
    print(f"The 'forbidden' pathway would lead to the {forbidden_pathway_product} isomer.")
    print("\n--- Conclusion ---")
    print("FMO theory predicts that the reaction proceeds exclusively through the allowed pathway.")
    print("Therefore, the predicted ratio of product A (cis) to product B (trans) is:")
    print(f"{ratio_A} : {ratio_B}")

# Run the analysis
predict_electrocyclic_product_ratio()