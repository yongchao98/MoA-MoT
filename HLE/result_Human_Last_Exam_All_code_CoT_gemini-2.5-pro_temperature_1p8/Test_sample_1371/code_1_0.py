import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to seat the conference members
    according to all the given constraints.
    """
    # Step 0: Define group composition based on the problem description.
    # Scientists: 12 total (6m, 6f). 2 male scientists are rowers.
    # -> Other Scientists: 4 male, 6 female.
    # Mathematicians: 4 total (2m, 2f). 1 male mathematician is a rower.
    # -> Other Mathematicians: 1 male, 2 female.
    # Classicists: 5 total. Cassie (f) + 4 others (2m, 2f).
    # Ethicists: 2 total (1m, 1f).
    # Rowers: 3 male total.

    s_other_male, s_other_female = 4, 6
    m_other_male, m_other_female = 1, 2
    
    # Pre-calculate factorials for efficiency
    fact = {n: math.factorial(n) for n in range(11)}

    # Step 1: Calculate ways for the internal arrangement of the SM-block.
    # The structure must be [10 other_S] - [3 Rowers] - [3 other_M].
    # N_internal = (ways to arrange 10 other_S with specific end) * 
    #              (ways to arrange 3 other_M with specific end) * 
    #              (ways to arrange 3 rowers).

    # Case 1 (Sf, Mf): End_S is female (6 choices), End_M is female (2 choices)
    ways_s_sf_end = s_other_female * fact[9]
    ways_m_mf_end = m_other_female * fact[2]
    n_int_1 = ways_s_sf_end * ways_m_mf_end * fact[3]

    # Case 2 (Sf, Mm): End_S is female (6 choices), End_M is male (1 choice)
    ways_m_mm_end = m_other_male * fact[2]
    n_int_2 = ways_s_sf_end * ways_m_mm_end * fact[3]

    # Case 3 (Sm, Mf): End_S is male (4 choices), End_M is female (2 choices)
    ways_s_sm_end = s_other_male * fact[9]
    n_int_3 = ways_s_sm_end * ways_m_mf_end * fact[3]

    # Case 4 (Sm, Mm): End_S is male (4 choices), End_M is male (1 choice)
    n_int_4 = ways_s_sm_end * ways_m_mm_end * fact[3]
    
    # Step 2: Calculate ways for the external arrangement of the 7 other people.
    # There are 7 seats around the SM-block: 2 adjacent and 5 non-adjacent ("middle").
    
    # N_ext for Case 1 (ends F, F): Both adjacent seats are OK for Cassie.
    # Place 4 other classicists in 5 middle seats: P(5,4).
    # Place Cassie + 2 Ethicists in the 3 remaining seats: 3!.
    p_5_4 = fact[5] // fact[1]
    n_ext_1 = p_5_4 * fact[3]

    # N_ext for Case 2 (ends F, M): One adjacent seat is OK for Cassie.
    # Place 4 other classicists in 5 middle seats: P(5,4).
    # Cassie has 2 choices (the OK adjacent seat + the leftover middle seat).
    # Place 2 Ethicists in the 2 remaining seats: 2!.
    n_ext_2 = p_5_4 * 2 * fact[2]

    # N_ext for Case 3 (ends M, F): Same as case 2.
    n_ext_3 = n_ext_2

    # N_ext for Case 4 (ends M, M): No adjacent seats are OK for Cassie.
    # All 5 classicists must go in the 5 middle seats: 5!.
    # Place 2 Ethicists in the 2 remaining (adjacent) seats: 2!.
    n_ext_4 = fact[5] * fact[2]

    # Step 3: Calculate the total ways for one block orientation (S-M) by summing the cases.
    w1 = n_int_1 * n_ext_1
    w2 = n_int_2 * n_ext_2
    w3 = n_int_3 * n_ext_3
    w4 = n_int_4 * n_ext_4
    w_sm_orientation = w1 + w2 + w3 + w4
    
    # Step 4: The final result must be multiplied by 2 for the two possible block orientations (S-M and M-S).
    total_ways = 2 * w_sm_orientation

    # Output the components of the final calculation as requested.
    print("The final calculation is 2 * (W_case1 + W_case2 + W_case3 + W_case4)\n")
    print(f"W_case1 (Female-Scientist-end, Female-Mathematician-end):")
    print(f"   Internal ways: {n_int_1}")
    print(f"   External ways: {n_ext_1}")
    print(f"   Subtotal: {n_int_1} * {n_ext_1} = {w1}\n")

    print(f"W_case2 (Female-Scientist-end, Male-Mathematician-end):")
    print(f"   Internal ways: {n_int_2}")
    print(f"   External ways: {n_ext_2}")
    print(f"   Subtotal: {n_int_2} * {n_ext_2} = {w2}\n")

    print(f"W_case3 (Male-Scientist-end, Female-Mathematician-end):")
    print(f"   Internal ways: {n_int_3}")
    print(f"   External ways: {n_ext_3}")
    print(f"   Subtotal: {n_int_3} * {n_ext_3} = {w3}\n")
    
    print(f"W_case4 (Male-Scientist-end, Male-Mathematician-end):")
    print(f"   Internal ways: {n_int_4}")
    print(f"   External ways: {n_ext_4}")
    print(f"   Subtotal: {n_int_4} * {n_ext_4} = {w4}\n")

    print(f"Total for one block orientation = {w1} + {w2} + {w3} + {w4} = {w_sm_orientation}\n")
    print(f"Total ways for both orientations = 2 * {w_sm_orientation} = {total_ways}")

solve_seating_arrangement()