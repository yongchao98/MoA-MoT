def calculate_bottcher_complexity():
    """
    This script calculates the Böttcher Molecular Complexity for the product
    of the Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Identify the product and its chemical properties.
    product_name = "cyclopentanecarboxylic acid"
    product_formula = "C6H10O2"
    print("Step 1: Identify the product of the reaction.")
    print(f"The Favorskii rearrangement of 2-chlorocyclohexanone results in a ring contraction,")
    print(f"producing {product_name}, with the chemical formula {product_formula}.\n")

    # Step 2: Calculate the parameters for the Böttcher Molecular Complexity formula.
    # The formula is: BMC = (Na * Nb) / Np
    print("Step 2: Calculate the parameters Na, Nb, and Np.")

    # Na: Total number of atoms from the formula C6H10O2
    num_c = 6
    num_h = 10
    num_o = 2
    Na = num_c + num_h + num_o
    print(f"Na (total atoms) = {num_c} (C) + {num_h} (H) + {num_o} (O) = {Na}")

    # Nb: Total number of bonds in the structure
    bonds_cc = 6        # 5 in ring, 1 to COOH
    bonds_ch = 9        # 9 on the ring
    bonds_co_double = 1 # C=O
    bonds_co_single = 1 # C-O
    bonds_oh = 1        # O-H
    Nb = bonds_cc + bonds_ch + bonds_co_double + bonds_co_single + bonds_oh
    print(f"Nb (total bonds) = {bonds_cc} (C-C) + {bonds_ch} (C-H) + {bonds_co_double} (C=O) + {bonds_co_single} (C-O) + {bonds_oh} (O-H) = {Nb}")

    # Np: Number of distinct elements (Carbon, Hydrogen, Oxygen)
    Np = 3
    print(f"Np (number of distinct elements: C, H, O) = {Np}\n")

    # Step 3: Compute the Böttcher Molecular Complexity (BMC).
    bmc_result = (Na * Nb) / Np

    print("Step 3: Calculate the Böttcher Molecular Complexity using the formula BMC = (Na * Nb) / Np.")
    print(f"Final Equation: BMC = ({Na} * {Nb}) / {Np}")
    print(f"The Böttcher Molecular Complexity is: {int(bmc_result)}")

calculate_bottcher_complexity()