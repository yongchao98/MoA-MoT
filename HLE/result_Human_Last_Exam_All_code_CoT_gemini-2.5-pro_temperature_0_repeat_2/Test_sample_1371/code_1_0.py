import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # Helper for factorial calculations
    def fact(n):
        return math.factorial(n)

    # --- Case 1: All 5 classicists sit together ---
    # Structure: [SM-Block] - E1 - [C-Group] - E2
    # The two ethicists can be placed in 2! ways around the SM-Block and C-Group.
    block_arrangements_case1 = fact(2)

    # Internal permutations of the SM-Block:
    # [10 S_other] - [2 S_rower] - [1 M_rower] - [3 M_other]
    # Arrange 10 other scientists: fact(10)
    # Arrange 2 rowing scientists: fact(2)
    # Arrange 3 other mathematicians: fact(3)
    internal_sm_block = fact(10) * fact(2) * fact(3)

    # Internal permutations of the 5 classicists
    internal_c_group = fact(5)

    # Total for Case 1
    ways_case_1 = block_arrangements_case1 * internal_sm_block * internal_c_group

    # --- Case 2: Cassie sits next to the SM-Block ---
    # This case is split into two sub-cases based on which end of the
    # SM-Block Cassie sits next to.

    # Sub-case 2a: Cassie sits next to the Scientist end.
    # This requires the scientist at that end to be female (6 choices).
    # Block arrangements: 2 ways (E1 or E2 next to SM-Block)
    block_arrangements_case2 = 2
    # Internal permutations of the 4 other classicists
    internal_c_reg_group = fact(4)

    # Conditional internal permutations of SM-Block (S-end must be female)
    # S-part: 6 (female choices for end) * fact(9) * fact(2) (rowers)
    # M-part: fact(3)
    internal_sm_cond_s = (6 * fact(9) * fact(2)) * fact(3)
    ways_case_2a = block_arrangements_case2 * internal_sm_cond_s * internal_c_reg_group

    # Sub-case 2b: Cassie sits next to the Mathematician end.
    # This requires the mathematician at that end to be female (2 choices).
    # Conditional internal permutations of SM-Block (M-end must be female)
    # S-part: fact(10) * fact(2)
    # M-part: 2 (female choices for end) * fact(2)
    internal_sm_cond_m = (fact(10) * fact(2)) * (2 * fact(2))
    ways_case_2b = block_arrangements_case2 * internal_sm_cond_m * internal_c_reg_group

    # Total ways is the sum of all cases
    total_ways = ways_case_1 + ways_case_2a + ways_case_2b

    # Print the final equation with all its components
    print("The total number of arrangements is the sum of three cases:")
    print("1. All classicists sit together.")
    print("2. Cassie sits by the scientist end of the main block.")
    print("3. Cassie sits by the mathematician end of the main block.")
    print("\nFinal Equation:")
    print(f"({block_arrangements_case1} * {internal_sm_block} * {internal_c_group}) + "
          f"({block_arrangements_case2} * {internal_sm_cond_s} * {internal_c_reg_group}) + "
          f"({block_arrangements_case2} * {internal_sm_cond_m} * {internal_c_reg_group})")
    print(f"= {ways_case_1} + {ways_case_2a} + {ways_case_2b}")
    print(f"= {total_ways}")

solve_seating_arrangement()
<<<13104424960>>>