import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the group at a circular table
    based on the given constraints.
    """

    # Part 1: Define group characteristics
    # There are 10 non-rowing-team scientists (6 female, 4 male)
    s_other_female = 6
    s_other_male = 4
    # There are 3 non-rowing-team mathematicians (2 female, 1 male)
    m_other_female = 2
    m_other_male = 1

    # Factorials needed for calculations
    # Factorial for arranging the 9 other scientists within the S-block
    f9 = math.factorial(9)
    # Factorial for arranging the 2 rowing scientists
    f2_s_rowers = math.factorial(2)
    # Factorial for arranging the 2 other mathematicians within the M-block
    f2_m_others = math.factorial(2)
    # Factorial for arranging the 5 people not seated next to the SM block
    f5 = math.factorial(5)

    # Part 2: Calculate internal SM-Block arrangements for each gender-end case

    # Number of ways to arrange the S-block part, ending with a specific gender
    # (Ways to choose end person) * (Ways to arrange other 9) * (Ways to arrange 2 rowers)
    s_arrangements_F_end = s_other_female * f9 * f2_s_rowers
    s_arrangements_M_end = s_other_male * f9 * f2_s_rowers

    # Number of ways to arrange the M-block part, ending with a specific gender
    # (Ways to choose end person) * (Ways to arrange other 2)
    m_arrangements_F_end = m_other_female * f2_m_others
    m_arrangements_M_end = m_other_male * f2_m_others

    # N_sm_internal for each case = (s_arrangements) * (m_arrangements)
    N_sm_internal_FF = s_arrangements_F_end * m_arrangements_F_end
    N_sm_internal_FM = s_arrangements_F_end * m_arrangements_M_end
    N_sm_internal_MF = s_arrangements_M_end * m_arrangements_F_end
    N_sm_internal_MM = s_arrangements_M_end * m_arrangements_M_end

    # Part 3: Calculate external arrangements for each case
    # There are 7 people to arrange around the SM block. The 2 adjacent seats
    # are special. The remaining 5 can be arranged in 5! ways.

    # Case 1 (Female, Female ends): Cassie and 2 Ethicists can sit adjacent.
    # Ways to choose 2 from 3 and arrange them: P(3,2)
    ext_ways_FF = math.perm(3, 2) * f5

    # Case 2 (Female, Male ends): Cassie can only sit on the Female side.
    # Adjacent seats are L (S-side) and R (M-side).
    # Ways to fill R (Male side) from {E1, E2}: 2
    # Ways to fill L (Female side) from {the other Ethicist, Cassie}: 2
    ext_ways_FM = (2 * 2) * f5

    # Case 3 (Male, Female ends): Symmetric to Case 2.
    ext_ways_MF = (2 * 2) * f5

    # Case 4 (Male, Male ends): Cassie cannot sit on either side.
    # Adjacent seats must be filled by the 2 Ethicists: P(2,2)
    ext_ways_MM = math.perm(2, 2) * f5

    # Part 4: Calculate total arrangements by summing the products for each case
    total_ways = (N_sm_internal_FF * ext_ways_FF +
                  N_sm_internal_FM * ext_ways_FM +
                  N_sm_internal_MF * ext_ways_MF +
                  N_sm_internal_MM * ext_ways_MM)

    # To show the final equation clearly, we find the constant factor.
    # Total = Factor * 9! * 5!
    # Factor = sum of [ (internal_coeffs) * (external_coeffs) ] for each case
    term1 = (s_other_female * f2_s_rowers * m_other_female * f2_m_others) * math.perm(3, 2)
    term2 = (s_other_female * f2_s_rowers * m_other_male * f2_m_others) * (2 * 2)
    term3 = (s_other_male * f2_s_rowers * m_other_female * f2_m_others) * (2 * 2)
    term4 = (s_other_male * f2_s_rowers * m_other_male * f2_m_others) * math.perm(2, 2)
    
    total_factor = term1 + term2 + term3 + term4
    
    print("The final calculation is broken down into a constant factor multiplied by the factorial arrangements of two groups of people.")
    print(f"The equation for the total number of arrangements is: {int(total_factor)} * {f9} * {f5}")
    print(f"Final Answer: {int(total_factor)} * {f9} * {f5} = {total_ways}")

solve_seating_arrangement()