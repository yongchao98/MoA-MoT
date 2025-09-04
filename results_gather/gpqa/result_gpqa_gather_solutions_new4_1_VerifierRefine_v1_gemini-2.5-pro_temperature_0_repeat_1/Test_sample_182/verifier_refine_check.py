def check_correctness():
    """
    Checks the correctness of the answer for the given chemistry problem.

    The question asks for the Index of Hydrogen Deficiency (IHD) of the product
    obtained when 2-formyl-5-vinylcyclohex-3-enecarboxylic acid is reacted
    with red phosphorus and excess of HI.
    """

    # --- Step 1: Analyze the starting material's structure ---
    # The IHD is the sum of rings and pi bonds.
    # 2-formyl-5-vinylcyclohex-3-enecarboxylic acid has:
    num_rings_start = 1          # from 'cyclohex-'
    num_pi_bonds_cc_start = 2    # from '-3-ene' and 'vinyl'
    num_pi_bonds_co_start = 2    # from 'formyl' and 'carboxylic acid'

    # --- Step 2: Analyze the reaction's effect ---
    # Red P + excess HI is a powerful reducing agent. It reduces all pi bonds
    # (both C=C and C=O) to single bonds but does not break the ring.

    # --- Step 3: Determine the product's structure ---
    # The product retains the ring but loses all pi bonds.
    num_rings_product = num_rings_start
    num_pi_bonds_product = 0

    # --- Step 4: Calculate the IHD of the product ---
    calculated_ihd_product = num_rings_product + num_pi_bonds_product

    # --- Step 5: Check against the provided answer ---
    # The options are A) 0, B) 5, C) 1, D) 3.
    # The provided answer is 'C', which corresponds to an IHD of 1.
    expected_ihd = 1

    if calculated_ihd_product == expected_ihd:
        return "Correct"
    else:
        reason = (
            f"The calculated Index of Hydrogen Deficiency (IHD) of the product is {calculated_ihd_product}, "
            f"but the provided answer 'C' corresponds to an IHD of {expected_ihd}.\n"
            "Reasoning:\n"
            "1. The starting material has 1 ring and 4 pi bonds (2 C=C, 2 C=O).\n"
            "2. The reaction with Red P + HI is a complete reduction, which saturates all pi bonds but preserves the ring.\n"
            f"3. Therefore, the product has {num_rings_product} ring and {num_pi_bonds_product} pi bonds.\n"
            f"4. The final IHD is the sum of rings and pi bonds, which is {num_rings_product} + {num_pi_bonds_product} = {calculated_ihd_product}."
        )
        return reason

# Run the check
result = check_correctness()
print(result)