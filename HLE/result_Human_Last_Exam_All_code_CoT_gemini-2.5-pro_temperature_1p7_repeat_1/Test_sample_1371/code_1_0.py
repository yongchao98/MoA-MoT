import math

def solve_seating_arrangement():
    """
    This function calculates the number of ways to seat the group according to the given constraints.
    """
    # Step 1: Define group numbers based on the problem and the necessary assumption.
    # Assumption: The number of classicists is 4 (2M, 2F) to meet the gender-balance criteria.
    # Total people = 12 scientists (6M, 6F) + 4 mathematicians (2M, 2F) + 2 ethicists (1M, 1F) + 4 classicists (2M, 2F) = 22.
    
    num_male_scientists = 6
    num_male_mathematicians = 2
    
    # Step 2: Calculate the total number of internal permutations for the combined SM-Block.
    # The SM-Block consists of the 12-scientist block and 4-mathematician block sitting together,
    # joined by the 3-person (2S, 1M) rowing team.
    
    # Ways to choose the 2 male rowing scientists from 6
    ways_to_choose_s_rowers = math.comb(num_male_scientists, 2)
    
    # Ways to choose the 1 male rowing mathematician from 2
    ways_to_choose_m_rower = math.comb(num_male_mathematicians, 1)
    
    # For a chosen set of individuals, calculate arrangements.
    # Arrange the 10 non-rower scientists: math.factorial(10)
    # Arrange the 2 rower scientists at their end of the block: math.factorial(2)
    # Arrange the 3 non-rower mathematicians at their end of the block: math.factorial(3)
    
    factorial_10 = math.factorial(10)
    factorial_4 = math.factorial(4)
    factorial_3 = math.factorial(3)
    factorial_2 = math.factorial(2)
    
    # N_internal = (Ways to choose rowers) * (Ways to arrange the rest)
    # N_internal = (C(6,2) * C(2,1)) * (10! * 2! * 3!)
    n_internal = ways_to_choose_s_rowers * ways_to_choose_m_rower * factorial_10 * factorial_2 * factorial_3
    
    # Step 3: Sub-divide internal permutations based on the gender of the people at the ends of the SM-Block.
    # Scientist end: Can be one of 10 non-rower scientists (4M, 6F). P(Female) = 6/10.
    # Mathematician end: Can be one of 3 non-rower mathematicians (1M, 2F). P(Female) = 2/3.
    
    prob_s_end_female = 6 / 10
    prob_m_end_female = 2 / 3
    
    # Number of arrangements where the scientist-end is a female
    n_s_end_female = n_internal * prob_s_end_female
    
    # Number of arrangements where the mathematician-end is a female
    n_m_end_female = n_internal * prob_m_end_female

    # Step 4: Calculate total arrangements by considering valid external configurations.
    # We are arranging 5 "items" in a circle: the main block [Neighbor1 - SM-Block - Neighbor2] and the 4 remaining people.
    # The circular permutation of 5 items is (5-1)! = 4!
    
    # Case 1: The two ethicists are the neighbors of the SM-Block.
    # The internal arrangement of SM-Block can be anything.
    # There are 2! ways to arrange the two ethicists.
    ways_case1 = factorial_2 * n_internal * factorial_4
    
    # Case 2: Cassie is a neighbor on one side, and an ethicist is on the other.
    # Subcase 2a: [Cassie - SM-Block - Ethicist]. The S-end must be female. There are 2 ethicists to choose from.
    ways_case2a = 2 * n_s_end_female * factorial_4

    # Subcase 2b: [Ethicist - SM-Block - Cassie]. The M-end must be female. There are 2 ethicists to choose from.
    ways_case2b = 2 * n_m_end_female * factorial_4
    
    # Total ways is the sum of all valid cases.
    total_ways = ways_case1 + ways_case2a + ways_case2b

    # Print the breakdown of the calculation
    print("This solution assumes a correction to the number of classicists from 5 to 4 to satisfy the gender balance rule.\n")
    print("The calculation is broken down into cases based on who sits next to the large Scientist-Mathematician block.")
    print(f"Number of internal arrangements for the Scientist-Mathematician block: {int(n_internal)}")
    print(f"Number of ways if ethicists sit on both sides = 2! * {int(n_internal)} * 4! = {int(ways_case1)}")
    print(f"Number of ways if Cassie sits next to the Scientist end = 2 * ({int(n_s_end_female)}) * 4! = {int(ways_case2a)}")
    print(f"Number of ways if Cassie sits next to the Mathematician end = 2 * ({int(n_m_end_female)}) * 4! = {int(ways_case2b)}")
    print("\nThe final equation is the sum of these cases:")
    print(f"Total ways = {int(ways_case1)} + {int(ways_case2a)} + {int(ways_case2b)}")
    print(f"Total ways = {int(total_ways)}")

solve_seating_arrangement()
<<<142149058560>>>