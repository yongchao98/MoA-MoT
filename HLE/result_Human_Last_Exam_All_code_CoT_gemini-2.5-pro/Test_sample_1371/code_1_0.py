import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a series of complex constraints.
    """
    # Step 1: Define constants based on the problem statement
    s_men = 6
    s_women = 6
    m_men = 2
    m_women = 2

    s_rowers_count = 2
    m_rowers_count = 1

    ethicists_count = 2
    
    # Calculate remaining people in Scientist and Mathematician groups after selecting rowers
    s_other_count = (s_men + s_women) - s_rowers_count
    s_other_women = s_women
    
    m_other_count = (m_men + m_women) - m_rowers_count
    m_other_women = m_women

    # Step 2: Calculate permutations for the SM super-block

    # Number of ways to choose the rowing team members from the available men
    n_choose_rowers = math.comb(s_men, s_rowers_count) * math.comb(m_men, m_rowers_count)

    # For a given set of rowers, calculate the permutations for the S-M and M-S block orders.
    # The arrangement is (other_S, rowers_S) linked to (rower_M, other_M).
    # Perms for [S][M] config: (10! * 2!) for S-block * (3!) for M-block.
    perms_sm_config = (math.factorial(s_other_count) * math.factorial(s_rowers_count)) * math.factorial(m_other_count)
    # Perms for [M][S] config is the same value.
    perms_ms_config = perms_sm_config

    # Total internal permutations is the sum for both orders, multiplied by ways to choose rowers.
    N_total_SM = n_choose_rowers * (perms_sm_config + perms_ms_config)

    # Step 3: Calculate conditional SM block permutations for Cassie's constraint

    # Fraction of arrangements with a female at the Scientist end
    frac_S_end_F = s_other_women / s_other_count
    # Fraction of arrangements with a female at the Mathematician end
    frac_M_end_F = m_other_women / m_other_count

    # N_SM_left_F: Permutations where the leftmost person of the super-block is female.
    # This occurs if order is [S][M] and S-end is female, OR order is [M][S] and M-end is female.
    N_SM_left_F = n_choose_rowers * (
        (perms_sm_config * frac_S_end_F) +
        (perms_ms_config * frac_M_end_F)
    )

    # N_SM_right_F: Permutations where the rightmost person is female.
    # This occurs if order is [S][M] and M-end is female, OR order is [M][S] and S-end is female.
    N_SM_right_F = n_choose_rowers * (
        (perms_sm_config * frac_M_end_F) +
        (perms_ms_config * frac_S_end_F)
    )

    # Step 4: Calculate external arrangement possibilities

    # There are 7 seats around the SM block: 2 adjacent (A) and 5 general (G).
    # The 4 other classicists must be placed in the 5 general seats.
    ways_place_C_prime = math.perm(5, 4)

    # The remaining 3 people (2 Ethicists, Cassie) go in the 3 remaining seats (1 G, 2 A).
    # Case A: Cassie is NOT adjacent (she takes the single remaining G seat).
    # The 2 Ethicists must take the 2 adjacent seats.
    ways_C_not_adj = ways_place_C_prime * 1 * math.factorial(ethicists_count)

    # Case B: Cassie IS adjacent. She can be on the Left or Right.
    # If Cassie takes an adjacent seat, the 2 Ethicists take the other two remaining seats.
    ways_C_adj = ways_place_C_prime * 1 * math.factorial(ethicists_count)

    # Step 5: Final Calculation
    term1 = ways_C_not_adj * N_total_SM
    term2 = ways_C_adj * N_SM_left_F
    term3 = ways_C_adj * N_SM_right_F

    total_ways = term1 + term2 + term3

    print("This problem is solved by considering the seating arrangements in several stages.")
    print("First, we treat the combined group of 12 scientists and 4 mathematicians as a single 'super-block'.")
    print("Second, we arrange this block and the 7 other individuals around the circular table.")
    print("Third, we apply all seating constraints, especially for the classicists.")
    print("\nThe final calculation is a sum of three parts, based on Cassie's position:")
    print("Total Ways = (Ways Cassie is NOT adjacent * Valid internal perms) + (Ways Cassie is on Left * Valid perms) + (Ways Cassie is on Right * Valid perms)")
    print("\nHere are the values for each part of the equation:")

    print(f"Ways for Cassie to NOT be adjacent to the super-block: {int(ways_C_not_adj)}")
    print(f"Total internal permutations for the super-block: {int(N_total_SM)}")
    print(f"Contribution from this case: {int(term1)}")

    print(f"\nWays for Cassie to be adjacent (e.g., on the Left): {int(ways_C_adj)}")
    print(f"Internal permutations with a Female on the Left end: {int(N_SM_left_F)}")
    print(f"Contribution from this case: {int(term2)}")

    print(f"\nWays for Cassie to be adjacent (e.g., on the Right): {int(ways_C_adj)}")
    print(f"Internal permutations with a Female on the Right end: {int(N_SM_right_F)}")
    print(f"Contribution from this case: {int(term3)}")

    print(f"\nFinal Equation:")
    print(f"{int(ways_C_not_adj)} * {int(N_total_SM)} + {int(ways_C_adj)} * {int(N_SM_left_F)} + {int(ways_C_adj)} * {int(N_SM_right_F)} = {int(total_ways)}")
    
    return int(total_ways)

final_answer = solve_seating_arrangement()
print(f"\n<<<Total number of ways to arrange the table: {final_answer}>>>")