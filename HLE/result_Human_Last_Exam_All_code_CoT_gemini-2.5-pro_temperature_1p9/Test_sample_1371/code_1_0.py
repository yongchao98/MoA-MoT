import math

def calculate_arrangements():
    """
    This function calculates the total number of seating arrangements based on the problem's constraints.
    """
    
    # Define a factorial function for convenience
    def factorial(n):
        return math.factorial(n)

    print("The problem is solved by considering two main configurations allowed by the constraints.")
    print("Let SM be the combined block of Scientists and Mathematicians.")
    print("Let C be the block of Classicists.")
    print("Let E1, E2 be the two Ethicists.")
    print("-" * 50)

    # --- Configuration A: Classicists separated by Ethicists ---
    print("Configuration A: The Classicists are fully separated from Scientists/Mathematicians.")
    print("The circular block arrangement must be SM - E1 - C - E2.")
    
    # 1. Arrangement of the blocks
    # We can swap the two distinct Ethicists (E1, E2), so there are 2 ways.
    n_blocks_A = 2
    print(f"Number of ways to arrange the blocks (swapping E1 and E2): {n_blocks_A}")

    # 2. Internal arrangement of the Classicist block (Block_C)
    # 5 classicists can be arranged in a line in 5! ways.
    internal_C_A = factorial(5)
    print(f"Internal arrangements for C block (5 people): 5! = {internal_C_A}")

    # 3. Internal arrangement of the Scientist-Mathematician block (Block_SM)
    # This block has an internal junction for the rowing team.
    # For Mathematicians (4 people): 1 male rower must be at the junction, others can be arranged.
    internal_M_A = factorial(4 - 1)
    # For Scientists (12 people): 2 male rowers must be at the junction, others can be arranged.
    internal_S_A = factorial(2) * factorial(12 - 2)
    internal_SM_A = internal_M_A * internal_S_A
    print(f"Internal arrangements for SM block:")
    print(f"  - M sub-block (1 rower at junction): (4-1)! = {internal_M_A}")
    print(f"  - S sub-block (2 rowers at junction): 2! * (12-2)! = {internal_S_A}")
    print(f"  - Total for SM block = {internal_M_A} * {internal_S_A} = {internal_SM_A}")
    
    ways_A = n_blocks_A * internal_C_A * internal_SM_A
    print(f"Total for Configuration A = {n_blocks_A} * {internal_C_A} * {internal_SM_A} = {ways_A}")
    print("-" * 50)
    
    # --- Configuration B: Cassie's Exception ---
    print("Configuration B: Classicists are adjacent to SM due to Cassie's exception.")
    print("The circular block arrangement must be SM - C - E1 - E2.")
    
    # 1. Arrangement of the blocks for this case
    # The composite block (SM-C), E1, and E2 can be arranged on a circle in (3-1)! = 2 ways.
    n_blocks_B = 2
    print(f"Number of ways to arrange the blocks (arranging SM-C, E1, E2): {n_blocks_B}")

    # 2. Internal arrangement of the Classicist block
    # Cassie must be at the junction, the other 4 can be arranged in 4! ways.
    internal_C_B = factorial(4)
    print(f"Internal arrangements for C block (Cassie at junction): (5-1)! = {internal_C_B}")
    
    # --- Sub-case B1: Junction at the Mathematician end of SM ---
    print("\nSub-case B1: The junction C-(M...S) requires a female Mathematician next to Cassie.")
    # M-block ways: 2 female mathematicians can be chosen. The rower is fixed at the other end.
    # The remaining 2 M's can be arranged in 2! ways.
    num_m_women = 2
    internal_M_B1 = num_m_women * factorial(4 - 2)
    # S-block ways are unchanged from Config A.
    internal_S_B1 = internal_S_A
    internal_SM_B1 = internal_M_B1 * internal_S_B1
    print(f"Internal arrangements for SM block in this sub-case:")
    print(f"  - M sub-block (female at junction): {num_m_women} * (4-2)! = {internal_M_B1}")
    print(f"  - S sub-block (unchanged): 2! * 10! = {internal_S_B1}")
    print(f"  - Total for SM block = {internal_M_B1} * {internal_S_B1} = {internal_SM_B1}")

    ways_B1 = n_blocks_B * internal_C_B * internal_SM_B1
    print(f"Total for Sub-case B1 = {n_blocks_B} * {internal_C_B} * {internal_SM_B1} = {ways_B1}")

    # --- Sub-case B2: Junction at the Scientist end of SM ---
    print("\nSub-case B2: The junction (M...S)-C requires a female Scientist next to Cassie.")
    # S-block ways: 6 female scientists can be chosen. Rowers are at the other end (2! ways).
    # Remaining 9 S's can be arranged in 9! ways.
    num_s_women = 6
    internal_S_B2 = num_s_women * factorial(12 - 1 - 2) * factorial(2)
    # M-block ways are unchanged from Config A.
    internal_M_B2 = internal_M_A
    internal_SM_B2 = internal_S_B2 * internal_M_B2
    print(f"Internal arrangements for SM block in this sub-case:")
    print(f"  - S sub-block (female at junction): {num_s_women} * (12-3)! * 2! = {internal_S_B2}")
    print(f"  - M sub-block (unchanged): 3! = {internal_M_B2}")
    print(f"  - Total for SM block = {internal_S_B2} * {internal_M_B2} = {internal_SM_B2}")

    ways_B2 = n_blocks_B * internal_C_B * internal_SM_B2
    print(f"Total for Sub-case B2 = {n_blocks_B} * {internal_C_B} * {internal_SM_B2} = {ways_B2}")
    print("-" * 50)
    
    # --- Final Calculation ---
    total_ways = ways_A + ways_B1 + ways_B2
    print("Total number of ways is the sum of all configurations:")
    print(f"Total = Ways(A) + Ways(B1) + Ways(B2)")
    print(f"Total = {ways_A} + {ways_B1} + {ways_B2} = {total_ways}")
    
    return total_ways

final_answer = calculate_arrangements()
print(f"\nFinal Answer: The total number of ways to arrange the table is {final_answer}.")
<<<13098516480>>>