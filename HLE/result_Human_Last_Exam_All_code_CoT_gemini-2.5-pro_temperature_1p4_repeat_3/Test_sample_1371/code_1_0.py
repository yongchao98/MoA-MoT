import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # Helper function for combinations
    def c(n, k):
        return math.comb(n, k)

    # Helper function for permutations
    def p(n, k):
        return math.perm(n, k)

    # Factorial function
    def f(n):
        return math.factorial(n)

    # --- Step 1: Define group characteristics ---
    s_men, s_women = 6, 6
    m_men, m_women = 2, 2
    num_ethicists = 2
    num_classicists = 5 # Cassie + 4 others
    s_rowers_m = 2
    m_rowers_m = 1

    # --- Step 2: Calculate internal arrangements for the Scientist block ---
    # The S-block is a line of 12 people. Two male rowers are at the 'connection' end.
    # We need to calculate arrangements based on the gender of the 'outer' end person.
    s_other_men = s_men - s_rowers_m
    s_other_women = s_women
    s_middle_positions = f(s_men + s_women - s_rowers_m - 1) # Factorial of 9 (12 - 2 rowers - 1 outer person)
    
    # Ways to choose and arrange the 2 male rowers from 6
    s_rower_config_ways = c(s_men, s_rowers_m) * f(s_rowers_m)
    
    # Arrangements with a FEMALE at the outer end
    ways_s_f_end = s_rower_config_ways * s_other_women * s_middle_positions
    # Arrangements with a MALE at the outer end
    ways_s_m_end = s_rower_config_ways * s_other_men * s_middle_positions
    ways_s_total = ways_s_f_end + ways_s_m_end

    # --- Step 3: Calculate internal arrangements for the Mathematician block ---
    # The M-block is a line of 4 people. One male rower is at the 'connection' end.
    m_other_men = m_men - m_rowers_m
    m_other_women = m_women
    m_middle_positions = f(m_men + m_women - m_rowers_m - 1) # Factorial of 2 (4 - 1 rower - 1 outer person)

    # Ways to choose the male rower from 2
    m_rower_config_ways = c(m_men, m_rowers_m)

    # Arrangements with a FEMALE at the outer end
    ways_m_f_end = m_rower_config_ways * m_other_women * m_middle_positions
    # Arrangements with a MALE at the outer end
    ways_m_m_end = m_rower_config_ways * m_other_men * m_middle_positions
    ways_m_total = ways_m_f_end + ways_m_m_end

    # --- Step 4: Calculate total arrangements by case ---

    # Case A: Neighbors of SM_super_block are the 2 Ethicists
    #   - Ways to place the 2 ethicists around the block: P(2,2)
    ways_place_ethicist_neighbors = p(num_ethicists, 2)
    #   - Ways to arrange the remaining 5 classicists: 5!
    ways_arrange_classicists = f(num_classicists)
    #   - Internal ways for the SM block (no end-person constraint)
    #     Product of total ways for S and M blocks, times 2 for S-M vs M-S orientation
    internal_ways_case_a = ways_s_total * ways_m_total * 2
    total_case_a = ways_place_ethicist_neighbors * ways_arrange_classicists * internal_ways_case_a

    # Case B: Neighbors are one Ethicist and Cassie
    #   - Ways to choose 1 of 2 ethicists and place them with Cassie: C(2,1)*P(2,2)
    ways_place_ethicist_cassie_neighbors = c(num_ethicists, 1) * p(2, 2)
    #   - Ways to arrange the other 5 people (1 E, 4 C): 5!
    ways_arrange_others_case_b = f(num_classicists)
    #   - Internal ways for the SM block, where one end (next to Cassie) MUST be female.
    #     We sum two sub-cases: (E-S...M-C) and (E-M...S-C)
    #     Sub-case E-S...M-C: S-block can be anything, M-block end must be female
    internal_ways_sm_cassie = ways_s_total * ways_m_f_end
    #     Sub-case E-M...S-C: M-block can be anything, S-block end must be female
    internal_ways_ms_cassie = ways_m_total * ways_s_f_end
    internal_ways_case_b = internal_ways_sm_cassie + internal_ways_ms_cassie
    total_case_b = ways_place_ethicist_cassie_neighbors * ways_arrange_others_case_b * internal_ways_case_b

    # --- Step 5: Final result and explanation ---
    total_arrangements = total_case_a + total_case_b

    print("The final number of arrangements is calculated by summing two main cases based on the neighbors of the combined Scientist-Mathematician block.")
    print("\n--- Case A: Ethicist-Ethicist Neighbors ---")
    print(f"Ways to place 2 Ethicists as neighbors = P(2,2) = {ways_place_ethicist_neighbors}")
    print(f"Ways to arrange the 5 Classicists = 5! = {ways_arrange_classicists}")
    print(f"Internal ways for S-M block = (Ways_S_Total * Ways_M_Total * 2) = ({ways_s_total} * {ways_m_total} * 2) = {internal_ways_case_a}")
    print(f"Subtotal for Case A = {ways_place_ethicist_neighbors} * {ways_arrange_classicists} * {internal_ways_case_a} = {total_case_a}")
    
    print("\n--- Case B: Ethicist-Cassie Neighbors ---")
    print(f"Ways to choose and place 1 Ethicist and Cassie as neighbors = C(2,1)*P(2,2) = {ways_place_ethicist_cassie_neighbors}")
    print(f"Ways to arrange the other 5 people = 5! = {ways_arrange_others_case_b}")
    print(f"Internal ways for S-M block (Cassie compatibility) = (Ways_S_Total * Ways_M_F_End) + (Ways_M_Total * Ways_S_F_End) = ({ways_s_total} * {ways_m_f_end}) + ({ways_m_total} * {ways_s_f_end}) = {internal_ways_case_b}")
    print(f"Subtotal for Case B = {ways_place_ethicist_cassie_neighbors} * {ways_arrange_others_case_b} * {internal_ways_case_b} = {total_case_b}")

    print("\n--- Final Calculation ---")
    print(f"Total Arrangements = Subtotal Case A + Subtotal Case B")
    print(f"Total Arrangements = {total_case_a} + {total_case_b} = {total_arrangements}")


solve_seating_arrangement()
<<<1421898240000>>>