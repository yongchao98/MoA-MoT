def solve_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Define the molecular formula of the starting material.
    # Starting Material: N-(p-methoxybenzyl)-5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one
    # C: 6(core) + 1(CF3) + 8(PMB) = 15
    # H: 5(core) + 9(PMB) = 14
    # F: 3(CF3) = 3
    # N: 1(core) = 1
    # O: 1(core) + 1(PMB) = 2
    # Formula: C15H14F3NO2
    initial_C = 15
    initial_H = 14
    initial_F = 3
    initial_N = 1
    initial_O = 2
    
    print("Starting Material Formula: C{}H{}F{}NO{}".format(initial_C, initial_H, initial_F, initial_O))
    print("---")

    # Step 2: Reaction 1 - PMB deprotection with CAN.
    # This removes the PMB group (C8H9O) and adds one hydrogen (H).
    # Net change: C: -8, H: -9 + 1 = -8, O: -1
    C_intermediate1 = initial_C - 8
    H_intermediate1 = initial_H - 8
    O_intermediate1 = initial_O - 1
    
    print("After Step 1 (PMB deprotection):")
    print(f"C = {initial_C} - 8 = {C_intermediate1}")
    print(f"H = {initial_H} - 9 + 1 = {H_intermediate1}")
    print(f"O = {initial_O} - 1 = {O_intermediate1}")
    print(f"Intermediate 1 Formula: C{C_intermediate1}H{H_intermediate1}F{initial_F}NO{O_intermediate1}")
    print("---")
    
    # Step 3: Reaction 2 - Hydrogenation of C=C bond with Pd/C, H2.
    # This adds one molecule of hydrogen (H2).
    # Net change: H: +2
    H_intermediate2 = H_intermediate1 + 2

    print("After Step 2 (Hydrogenation):")
    print(f"H = {H_intermediate1} + 2 = {H_intermediate2}")
    print(f"Intermediate 2 Formula: C{C_intermediate1}H{H_intermediate2}F{initial_F}NO{O_intermediate1}")
    print("---")
    
    # Step 4: Reaction 3 - Hydrolysis of the lactam with acid.
    # This adds one molecule of water (H2O).
    # Net change: H: +2, O: +1
    product_C = C_intermediate1
    product_H = H_intermediate2 + 2
    product_F = initial_F
    product_N = initial_N
    product_O = O_intermediate1 + 1
    
    print("After Step 3 (Lactam Hydrolysis):")
    print(f"H = {H_intermediate2} + 2 = {product_H}")
    print(f"O = {O_intermediate1} + 1 = {product_O}")
    print("---")
    
    # Final Result
    print("The molecular formula of the final product is:")
    print(f"C{product_C}H{product_H}F{product_F}NO{product_O}")

solve_molecular_formula()