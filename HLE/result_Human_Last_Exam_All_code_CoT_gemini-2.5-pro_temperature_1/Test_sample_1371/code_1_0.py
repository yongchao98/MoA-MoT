import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the people at the conference table
    based on the given constraints.
    """
    print("This problem is solved by breaking it into two parts: the external arrangement of groups and the internal arrangement of the main block.\n")

    # Part 1: External Arrangement Calculation
    print("--- Part 1: External Arrangement ---")
    print("We are arranging 8 'items' in a circle: a Scientist-Mathematician block (Block_SM), 2 Ethicists, 4 regular Classicists, and Cassie.")
    print("The 4 regular Classicists cannot sit next to Block_SM. The 2 Ethicists and Cassie can.")
    
    # The two seats next to the Block_SM must be filled by 2 of the 3 people who can sit there.
    # This is the number of permutations of choosing 2 from 3.
    p_3_2 = math.perm(3, 2)
    print(f"Number of ways to place neighbors for Block_SM (P(3, 2)): {p_3_2}")

    # The remaining 5 people/items can be arranged in the remaining 5 seats.
    fact_5 = math.factorial(5)
    print(f"Number of ways to arrange the remaining 5 people (5!): {fact_5}")

    # Total external arrangements
    n_external = p_3_2 * fact_5
    print(f"Total External Arrangements = {p_3_2} * {fact_5} = {n_external}\n")

    # Part 2: Internal Arrangement of the Scientist-Mathematician Block
    print("--- Part 2: Internal Arrangement of Block_SM ---")
    print("This block contains 12 scientists (6M, 6W) and 4 mathematicians (2M, 2W).")
    print("A rowing team of 2 male scientists and 1 male mathematician must sit together at the junction of the two groups.")

    # Ways to choose the rowing team members
    c_2_1 = math.comb(2, 1) # Choose 1 of 2 male mathematicians
    c_6_2 = math.comb(6, 2) # Choose 2 of 6 male scientists
    ways_to_choose_rowers = c_2_1 * c_6_2
    print(f"Ways to choose 1 male mathematician from 2 (C(2, 1)): {c_2_1}")
    print(f"Ways to choose 2 male scientists from 6 (C(6, 2)): {c_6_2}")
    print(f"Total ways to choose the rowing team = {c_2_1} * {c_6_2} = {ways_to_choose_rowers}")

    # Ways to arrange the members within the block
    fact_3 = math.factorial(3) # Arrange the other 3 mathematicians
    fact_2 = math.factorial(2) # Arrange the 2 scientist rowers among themselves
    fact_10 = math.factorial(10) # Arrange the other 10 scientists
    print(f"Ways to arrange the other 3 mathematicians (3!): {fact_3}")
    print(f"Ways to arrange the 2 scientist rowers (2!): {fact_2}")
    print(f"Ways to arrange the other 10 scientists (10!): {fact_10}")
    
    ways_to_arrange_block_members = fact_3 * fact_2 * fact_10
    
    # Total internal arrangements
    n_internal = ways_to_choose_rowers * ways_to_arrange_block_members
    print(f"Total Internal Arrangements = {ways_to_choose_rowers} * ({fact_3} * {fact_2} * {fact_10}) = {n_internal}\n")

    # Final Calculation
    print("--- Final Calculation ---")
    total_arrangements = n_external * n_internal
    
    print("Total Arrangements = (External Arrangements) * (Internal Arrangements)")
    print(f"Total = ({p_3_2} * {fact_5}) * ({c_2_1} * {c_6_2} * {fact_3} * {fact_2} * {fact_10})")
    print(f"Total = {n_external} * {n_internal}")
    print(f"The total number of ways to arrange the table is: {total_arrangements}")
    
    return total_arrangements

final_answer = solve_seating_arrangement()
print(f"\n<<<{final_answer}>>>")
