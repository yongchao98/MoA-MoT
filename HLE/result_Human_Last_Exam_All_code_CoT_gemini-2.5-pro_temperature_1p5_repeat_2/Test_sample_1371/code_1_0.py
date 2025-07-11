import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # Step 1: Define constants from the problem description
    num_scientists = 12
    num_mathematicians = 4
    num_classicists = 5
    num_ethicists = 2

    # Gender distribution (equal for each profession, Cassie's gender is unknown but irrelevant for the counts)
    # Other scientists (total 12 - 2 rowers = 10): 4 male, 6 female
    # Other mathematicians (total 4 - 1 rower = 3): 1 male, 2 female
    s_other_female = 6
    s_other_total = 10
    m_other_female = 2
    m_other_total = 3

    # Step 2: Calculate the internal permutations of the Scientist-Mathematician (SM) Block.
    # This block has 16 people.
    # The two groups (S and M) must be adjacent. This gives a factor of 2 (S-M or M-S).
    # The 2 scientist rowers can be arranged in 2! ways at the seam.
    # The remaining 10 scientists can be arranged in 10! ways.
    # The remaining 3 mathematicians can be arranged in 3! ways.
    
    f_10 = math.factorial(10)
    f_3 = math.factorial(3)
    f_2 = math.factorial(2)
    f_5 = math.factorial(5)

    W_SM = 2 * f_2 * f_10 * f_3

    # Step 3: Calculate the number of arrangements for each scenario.

    # Scenario A: Both neighbors of the SM-Block are Ethicists (E-SM-E).
    # Treat (E-SM-E) as one block. The 2 ethicists can be on either side (2! ways).
    # Arrange this E-SM-E block and the 5 classicists in a circle: (1 + 5 - 1)! = 5! ways.
    # Ways_A = (permutations of E-SM-E block) * (circular permutations of items)
    perms_EBE_block = f_2 * W_SM
    ways_A = f_5 * perms_EBE_block

    # Scenario B: One neighbor is an Ethicist, the other is Cassie.
    # This requires calculating the number of SM-Block permutations where the person
    # next to Cassie is female.

    # Probability that the person at the Scientist-end of the block is female.
    prob_s_end_female = s_other_female / s_other_total
    # Probability that the person at the Mathematician-end of the block is female.
    prob_m_end_female = m_other_female / m_other_total
    
    # Valid permutations for Cassie sitting next to the S-side.
    W_SM_S_F = W_SM * prob_s_end_female
    # Valid permutations for Cassie sitting next to the M-side.
    W_SM_M_F = W_SM * prob_m_end_female

    # There are 2 choices for which ethicist acts as the buffer.
    # The remaining items to arrange are the (Cassie-SM-E) block, the other ethicist,
    # and the 4 other classicists. In a circle: (1 + 1 + 4 - 1)! = 5! ways.

    ways_B1_Cassie_at_M_side = 2 * f_5 * W_SM_M_F
    ways_B2_Cassie_at_S_side = 2 * f_5 * W_SM_S_F

    # Step 4: Sum the results for the total number of arrangements.
    total_ways = ways_A + ways_B1_Cassie_at_M_side + ways_B2_Cassie_at_S_side
    
    # --- Final Output ---
    print("This problem is solved by partitioning it into three scenarios based on who sits next to the combined Scientist-Mathematician (SM) block.")
    print("\n1. Internal permutations of the SM-Block (W_SM):")
    print(f"   W_SM = 2 * 2! * 10! * 3! = {W_SM:,.0f}")
    
    print("\n2. Calculating permutations for each scenario:")
    print("\n   Scenario A: Both neighbors are Ethicists (E-SM-E).")
    print(f"   Ways_A = 5! * (2! * W_SM) = {f_5} * ({f_2} * {W_SM:,.0f}) = {ways_A:,.0f}")
    
    print("\n   Scenario B1: Cassie sits next to the Mathematician side (E-SM-C).")
    print(f"   Valid W_SM for M-side female = W_SM * (2/3) = {W_SM_M_F:,.0f}")
    print(f"   Ways_B1 = 2 * 5! * (Valid W_SM) = 2 * {f_5} * {W_SM_M_F:,.0f} = {ways_B1_Cassie_at_M_side:,.0f}")

    print("\n   Scenario B2: Cassie sits next to the Scientist side (C-SM-E).")
    print(f"   Valid W_SM for S-side female = W_SM * (6/10) = {W_SM_S_F:,.0f}")
    print(f"   Ways_B2 = 2 * 5! * (Valid W_SM) = 2 * {f_5} * {W_SM_S_F:,.0f} = {ways_B2_Cassie_at_S_side:,.0f}")

    print("\n3. Total number of arrangements is the sum of all scenarios:")
    print(f"   Total Ways = Ways_A + Ways_B1 + Ways_B2")
    print(f"   Total Ways = {ways_A:,.0f} + {ways_B1_Cassie_at_M_side:,.0f} + {ways_B2_Cassie_at_S_side:,.0f}")
    print(f"   Total Ways = {total_ways:,.0f}")

solve_seating_arrangement()