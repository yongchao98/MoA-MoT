import math

def solve_conference_seating():
    """
    This function calculates the number of ways to arrange the members at the circular table
    based on the given constraints.
    """

    # Step 1: Define the constants from the problem description.
    num_scientists = 12
    num_mathematicians = 4
    num_ethicists = 2
    num_classicists = 5
    
    # Gender distribution
    female_scientists = 6
    female_mathematicians = 2
    
    # The problem is solved by fixing the large 16-person Scientist-Mathematician (SM) block
    # and then arranging the remaining 7 people (2 Ethicists, 5 Classicists) around it.
    # The 4 regular classicists (RCs) cannot sit next to the SM block.
    # This leaves the 2 seats next to the SM block to be filled by the 2 Ethicists and/or Cassie.

    # Let's calculate the number of arrangements for the two main cases for the SM block's neighbors.

    # --- Case 1: The two neighbors of the SM block are the two Ethicists. ---
    
    # Ways to arrange the two Ethicists in the two neighbor seats:
    ways_place_ethicists = math.factorial(num_ethicists) # 2!
    
    # The remaining 5 people are all classicists, arranged in the 5-seat gap:
    ways_arrange_classicists_gap = math.factorial(num_classicists) # 5!
    
    # Internal permutations of the SM block (unconstrained by Cassie):
    # S-block: 10 non-rowers (10!) and 2 rowers at the end (2!).
    # M-block: 3 non-rowers (3!) and 1 rower at the end (1!).
    internal_sm_unconstrained = (math.factorial(10) * math.factorial(2)) * math.factorial(3)
    
    # Total for Case 1:
    total_case_1 = ways_place_ethicists * ways_arrange_classicists_gap * internal_sm_unconstrained
    
    # --- Case 2: The two neighbors are one Ethicist and Cassie. ---
    
    # Ways to choose which of the 2 Ethicists is a neighbor:
    ways_choose_ethicist = 2
    
    # Ways to place Cassie and the chosen Ethicist in the 2 neighbor seats:
    ways_place_cassie_and_ethicist = math.factorial(2) # 2!
    
    # The remaining 5 people (4 RCs and the other Ethicist) are arranged in the 5-seat gap:
    ways_arrange_others_gap = math.factorial(5) # 5!
    
    # Internal permutations of the SM block are now constrained by Cassie.
    
    # Subcase 2a: Cassie is at the Scientist end.
    # The scientist at her end must be one of the 6 females.
    # S-block: 6 choices for female at end * 9! for the rest * 2! for rowers at other end.
    # M-block: 3! for the rest.
    internal_sm_cassie_at_s_end = (female_scientists * math.factorial(9) * math.factorial(2)) * math.factorial(3)

    # Subcase 2b: Cassie is at the Mathematician end.
    # The mathematician at her end must be one of the 2 females.
    # S-block: 10! * 2! for the rest.
    # M-block: 2 choices for female at end * 2! for the rest * 1 for rower at other end.
    internal_sm_cassie_at_m_end = (math.factorial(10) * math.factorial(2)) * (female_mathematicians * math.factorial(2))

    # Total for Case 2:
    # We have 2 choices for the Ethicist. For each choice, Cassie can be at the S-end or M-end.
    # The arrangement of neighbors can be (C, E) or (E, C).
    # Total = (ways to choose E) * (ways to arrange others in gap) * [(perms for C at S-end) + (perms for C at M-end)]
    # The factor of 2 for (C,E) vs (E,C) is implicitly handled by summing the S-end and M-end cases.
    total_case_2 = ways_choose_ethicist * ways_arrange_others_gap * (internal_sm_cassie_at_s_end + internal_sm_cassie_at_m_end)

    # --- Final Calculation ---
    # The total number of arrangements is the sum of Case 1 and Case 2.
    # Total = [2! * 5! * (10!*2!*3!)] + [2 * 5! * ((6*9!*2!*3!) + (10!*2!*4))]
    # This simplifies to 240 * 272 * 9!
    
    total_ways = total_case_1 + total_case_2
    
    # To show the simplified calculation as requested:
    # 240 comes from 2 * 5! (from Case 1) or 2 * 5! (from Case 2)
    # 272 * 9! comes from the sum of the internal permutation parts.
    # [ (10!*2!*3!) + ((6*9!*2!*3!) + (10!*2!*4)) ] / (10!*12) ... it's complex.
    # Let's use the simplified factors derived during the thinking process.
    
    factor1 = 240
    factor2 = 272
    factor3 = math.factorial(9)
    
    calculated_total = factor1 * factor2 * factor3
    
    print("The problem can be broken down into cases based on who sits next to the combined Scientist-Mathematician block.")
    print("Summing the possibilities for these cases leads to a simplified final equation.")
    print(f"The final calculation is {factor1} * {factor2} * {factor3} = {calculated_total}")

solve_conference_seating()
<<<23688806400>>>