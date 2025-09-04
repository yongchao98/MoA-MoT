def check_chemistry_answer():
    """
    Checks the correctness of the answer for the IHD of the product of a given reaction.

    Question: What is the index of hydrogen deficiency of the product obtained when
    2-formyl-5-vinylcyclohex-3-enecarboxylic acid is reacted with red phosphorus and excess of HI?

    Answer: D) 1
    """

    # Step 1: Define the sources of IHD in the reactant molecule.
    # Reactant: 2-formyl-5-vinylcyclohex-3-enecarboxylic acid
    reactant_ihd_sources = {
        "ring": 1,
        "c=c_in_ring": 1,
        "c=c_in_vinyl_group": 1,
        "c=o_in_formyl_group": 1,
        "c=o_in_carboxylic_acid": 1
    }
    ihd_reactant = sum(reactant_ihd_sources.values())

    # Step 2: Define the effect of the reaction (Red P + excess HI).
    # This is a strong reduction that saturates pi bonds but does not open stable rings.
    # We identify which sources of IHD are removed by the reaction.
    pi_bonds_reduced = (
        reactant_ihd_sources["c=c_in_ring"] +
        reactant_ihd_sources["c=c_in_vinyl_group"] +
        reactant_ihd_sources["c=o_in_formyl_group"] +
        reactant_ihd_sources["c=o_in_carboxylic_acid"]
    )
    rings_broken = 0 # Red P/HI does not open stable cycloalkane rings.

    # Step 3: Calculate the IHD of the final product.
    ihd_product = ihd_reactant - pi_bonds_reduced - rings_broken

    # Step 4: Check if the calculated IHD matches the provided answer.
    # The provided answer is D, which corresponds to an IHD of 1.
    correct_answer_ihd = 1

    if ihd_product == correct_answer_ihd:
        return "Correct"
    else:
        error_message = (
            f"The answer is incorrect.\n"
            f"Reasoning:\n"
            f"1. The IHD of the reactant (2-formyl-5-vinylcyclohex-3-enecarboxylic acid) is calculated by summing its unsaturation elements: 1 (ring) + 1 (C=C in ring) + 1 (C=C in vinyl) + 1 (C=O in formyl) + 1 (C=O in acid) = {ihd_reactant}.\n"
            f"2. The reaction with Red Phosphorus and excess HI is a powerful reduction. It reduces all 4 pi bonds (two C=C and two C=O) to saturated single bonds.\n"
            f"3. The stable cyclohexane ring is not opened by this reaction.\n"
            f"4. Therefore, the IHD of the product is the IHD of the reactant minus the IHD from the reduced pi bonds: {ihd_reactant} - {pi_bonds_reduced} = {ihd_product}.\n"
            f"The calculated IHD is {ihd_product}, but the answer given corresponds to an IHD of {correct_answer_ihd}."
        )
        return error_message

# Execute the check
result = check_chemistry_answer()
print(result)