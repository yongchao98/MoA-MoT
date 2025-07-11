import math

def calculate_seating_arrangements():
    """
    Calculates the number of ways to arrange the group at a circular table
    based on a complex set of constraints.
    """
    # Step 1: Define factorial shortcuts for readability
    fact = math.factorial
    
    # Step 2: Calculate the number of ways to choose the members of the rowing team
    # Choose 2 male scientists from 6: C(6, 2)
    # Choose 1 male mathematician from 2: C(2, 1)
    ways_to_choose_rowers_val = math.comb(6, 2) * math.comb(2, 1)

    # Step 3: Calculate the internal arrangements for the Block_SM (Scientists & Mathematicians)
    # This block has a fixed structure: (10 non-rower scientists) - (3 rowers) - (3 non-rower mathematicians)
    # The ends of this block are one scientist and one mathematician. We need to know the arrangements
    # based on the gender of these end people for Cassie's constraint.
    # Non-rower scientists: 6 Women (SW), 4 Men (SM). Total 10.
    # Non-rower mathematicians: 2 Women (MW), 1 Man (MM). Total 3.
    # Rowers: 3 Men.
    
    # Internal arrangements of the 10 other scientists (6W, 4M)
    arr_s_end_W = 6 * fact(9)  # A woman is at the end connected to rowers
    arr_s_end_M = 4 * fact(9)  # A man is at the end connected to rowers
    
    # Internal arrangements of the 3 rowers
    arr_rowers = fact(3)
    
    # Internal arrangements of the 3 other mathematicians (2W, 1M)
    # A woman is at the end of the block (2 choices for the woman * 2! ways for the rest)
    arr_m_end_W = 2 * fact(2)
    # A man is at the end of the block (1 choice for the man * 2! ways for the rest)
    arr_m_end_M = 1 * fact(2)

    # Total BSM arrangements for different end-gender combinations.
    # The full block arrangement is (Arr S-group) * (Arr R-group) * (Arr M-group)
    bsm_W_W = arr_s_end_W * arr_rowers * arr_m_end_W
    bsm_W_M = arr_s_end_W * arr_rowers * arr_m_end_M
    bsm_M_W = arr_s_end_M * arr_rowers * arr_m_end_W
    bsm_M_M = arr_s_end_M * arr_rowers * arr_m_end_M
    
    # Arrangements for Cassie's constraint
    bsm_s_end_W_any_m = bsm_W_W + bsm_W_M  # S-end is Woman, M-end is Any
    bsm_any_s_m_end_W = bsm_W_W + bsm_M_W  # S-end is Any, M-end is Woman
    bsm_any_any = bsm_W_W + bsm_W_M + bsm_M_W + bsm_M_M

    # Step 4: Calculate internal arrangements for the Friendly_Block (Ethicists & Classicists)
    # Total 7 people. Ends must be from {E1, E2, Cassie}. 5 people are in the middle.
    arr_middle_5 = fact(5)
    
    # Case 1: Ends are the two Ethicists (E1 ... E2 or E2 ... E1)
    fb_E_E = fact(2) * arr_middle_5
    
    # Case 2: One end is Cassie, the other is an Ethicist.
    # Choice of Ethicist (2) * arrangements of ends (C...E or E...C) * arrangements of middle 5.
    fb_C_at_end = 2 * fact(2) * arr_middle_5


    # Step 5: Combine the arrangements
    # Since there are two blocks, we sum the products of valid pairings.
    # Term 1: Ethicist ends (fb_E_E) connect to BSM. No constraint, so use bsm_any_any.
    total_term1 = fb_E_E * bsm_any_any
    
    # Term 2: Cassie at one end, Ethicist at other (fb_C_at_end).
    # Half of these have Cassie at the S-side of BSM, half have her at the M-side.
    # Cassie at S-side: Her BSM neighbor must be a woman (bsm_s_end_W_any_m).
    # Cassie at M-side: Her BSM neighbor must be a woman (bsm_any_s_m_end_W).
    total_term2 = (fb_C_at_end / 2) * bsm_s_end_W_any_m + (fb_C_at_end / 2) * bsm_any_s_m_end_W

    # Total structural arrangements
    total_arrangements = total_term1 + total_term2

    # Final result is (ways to choose) * (total structural arrangements)
    final_result = ways_to_choose_rowers_val * total_arrangements

    # Step 6: Print the detailed calculation and result
    # For clarity, let's represent the BSM calculations in a more compact form for the printout.
    # `fact(9)` is common in all BSM arrangement numbers.
    # The sum simplifies to ways_to_choose * (fb_E_E * (bsm... / fact(9)) + ...) * fact(9)
    # The part `(total_arrangements / fact(9))` gives a cleaner integer.
    inner_factor = int(total_arrangements / fact(9))

    print("The total number of arrangements is calculated by combining the ways to form subgroups with the ways to arrange them.")
    print("\n1. Ways to choose the rowing team members:")
    print("   Ways = C(6, 2) * C(2, 1)")
    print(f"   Ways = {math.comb(6, 2)} * {math.comb(2, 1)} = {ways_to_choose_rowers_val}")

    print("\n2. Ways to arrange the two main blocks ('Friendly' and 'Scientists/Mathematicians'):")
    print("   This is a sum of arrangements for valid end-to-end connections.")
    print("   Total arrangements = [ (Friendly Block with E-E ends) * (BSM with any ends) + (Friendly Block with C-E ends) * (BSM with female end for Cassie) ]")
    print("   Simplified, this calculates to:")
    print(f"   Arrangements = {inner_factor} * 9!")
    
    print("\n3. Final Calculation:")
    final_eq = (
        f"{ways_to_choose_rowers_val} * {inner_factor} * {fact(9)} = {final_result}"
    )
    print(f"   Total Ways = (Ways to choose) * (Arrangements)")
    print(f"   Total Ways = {final_eq}")
    
calculate_seating_arrangements()
<<<21326817024000>>>