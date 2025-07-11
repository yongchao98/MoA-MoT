import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the guests at the circular table
    based on the given constraints.
    """

    # Step 1: Define group numbers from the problem description
    # Scientists: 12 total -> 2 male rowers, 4 male non-rowers, 6 female non-rowers
    n_total_scientists = 12
    n_female_scientists = 6
    n_non_rower_scientists = 10 # 12 total - 2 rowers

    # Mathematicians: 4 total -> 1 male rower, 1 male non-rower, 2 female non-rowers
    n_total_mathematicians = 4
    n_female_mathematicians = 2
    n_non_rower_mathematicians = 3 # 4 total - 1 rower

    # Other groups
    n_rowers = 3
    n_ethicists = 2
    n_classicists_prime = 4 # Classicists excluding Cassie

    # Step 2: Calculate permutations for the Super-Block SM, considering Cassie's preference.
    # Cassie must sit next to a female scientist or a female mathematician at the ends of the SM block.

    # Case 1: Cassie sits next to the Scientist end. The person at that end must be a female scientist.
    # Ways = (choices for female scientist at the end) * (permutations of remaining scientists) *
    #        (permutations of rowers) * (permutations of mathematicians)
    ways_sm_s_end = n_female_scientists * math.factorial(n_non_rower_scientists - 1) * math.factorial(n_rowers) * math.factorial(n_non_rower_mathematicians)

    # Case 2: Cassie sits next to the Mathematician end. The person at that end must be a female mathematician.
    # Ways = (choices for female math at the end) * (permutations of remaining mathematicians) *
    #        (permutations of rowers) * (permutations of scientists)
    ways_sm_m_end = n_female_mathematicians * math.factorial(n_non_rower_mathematicians - 1) * math.factorial(n_rowers) * math.factorial(n_non_rower_scientists)

    # Total permutations for the SM block where Cassie's preference is met.
    total_ways_sm = ways_sm_s_end + ways_sm_m_end

    # Step 3: Calculate permutations for the other blocks.
    ways_ethicists = math.factorial(n_ethicists)
    ways_classicists_prime = math.factorial(n_classicists_prime)

    # Step 4: Calculate total arrangements.
    # There are 2 ways to arrange the blocks (SM-E-C'-Cassie and SM-Cassie-C'-E).
    num_block_arrangements = 2
    
    total_arrangements = num_block_arrangements * total_ways_sm * ways_ethicists * ways_classicists_prime
    
    # Step 5: Print the breakdown and the final answer.
    print("The final calculation is broken down as follows:")
    print("Total Ways = (Block Arrangements) * (SM Permutations) * (Ethicist Permutations) * (Classicist Permutations)")
    print("\nLet's break down the SM Permutations (the most complex part):")
    print("SM_Perms = (Cassie next to female Scientist) + (Cassie next to female Mathematician)")
    print(f"  - Part 1 (Cassie next to Scientist): {n_female_scientists} * {n_non_rower_scientists - 1}! * {n_rowers}! * {n_non_rower_mathematicians}! = {ways_sm_s_end}")
    print(f"  - Part 2 (Cassie next to Mathematician): {n_female_mathematicians} * {n_non_rower_mathematicians - 1}! * {n_rowers}! * {n_non_rower_scientists}! = {ways_sm_m_end}")
    print(f"Total SM Permutations = {ways_sm_s_end} + {ways_sm_m_end} = {total_ways_sm}")
    
    print("\nPutting it all together:")
    final_equation = f"{num_block_arrangements} * ({ways_sm_s_end} + {ways_sm_m_end}) * {n_ethicists}! * {n_classicists_prime}!"
    print(f"Total Ways = {final_equation}")
    print(f"Total Ways = {num_block_arrangements} * {total_ways_sm} * {ways_ethicists} * {ways_classicists_prime}")
    print(f"\nFinal calculated number of ways: {total_arrangements}")


solve_seating_arrangement()
<<<15885434880>>>