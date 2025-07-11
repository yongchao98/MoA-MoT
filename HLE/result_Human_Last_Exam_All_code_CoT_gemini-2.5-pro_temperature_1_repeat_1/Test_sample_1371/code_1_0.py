import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to arrange the members at a circular table
    based on a set of complex constraints.
    """

    # PART 1: Internal Arrangements of the Scientist-Mathematician Super-Block
    # The block consists of 12 scientists and 4 mathematicians.
    # A rowing team (2 male scientists, 1 male math) forces these two groups to be adjacent.
    # The block has a Scientist end and a Mathematician end.
    #
    # Scientists: 12 total (6M, 6W). The 2 rowing scientists are male.
    # Other scientists: 10 people (4M, 6W).
    # Mathematicians: 4 total (2M, 2W). The 1 rowing mathematician is male.
    # Other mathematicians: 3 people (1M, 2W).

    # Number of ways to arrange the 'other' 9 scientists behind the Scientist at the end.
    fact_9 = math.factorial(9)
    # Number of ways to arrange the 2 rowing scientists at their end of the scientist line.
    fact_2 = math.factorial(2)
    # Number of ways to arrange the 2 'middle' mathematicians.
    fact_2_math = math.factorial(2)

    # Calculate internal arrangements based on gender of the two endpoints.
    # Let S_end be the scientist at one end and M_end be the mathematician at the other.

    # N_FF: S_end is Female, M_end is Female
    # Ways to choose a female S_end: 6
    # Ways to choose a female M_end: 2
    ways_s_end_F = 6 * fact_9 * fact_2
    ways_m_end_F = 2 * fact_2_math
    N_FF = ways_s_end_F * ways_m_end_F

    # N_FM: S_end is Female, M_end is Male
    # Ways to choose a male M_end: 1
    ways_m_end_M = 1 * fact_2_math
    N_FM = ways_s_end_F * ways_m_end_M

    # N_MF: S_end is Male, M_end is Female
    # Ways to choose a male S_end from non-rowers: 4
    ways_s_end_M = 4 * fact_9 * fact_2
    N_MF = ways_s_end_M * ways_m_end_F

    # N_MM: S_end is Male, M_end is Male
    N_MM = ways_s_end_M * ways_m_end_M

    # Total internal arrangements if there are no endpoint constraints.
    N_total = N_FF + N_FM + N_MF + N_MM

    # Internal arrangements where the Scientist-end must be female (for Cassie).
    N_S_female = N_FF + N_FM
    # Internal arrangements where the Mathematician-end must be female (for Cassie).
    N_M_female = N_FF + N_MF

    # PART 2: External Arrangements around the Circular Table
    # The entities are: the Super-Block (U_SM), 2 Ethicists, and 5 Classicists.
    # The 4 'regular' classicists cannot sit next to the Super-Block.
    # The neighbors must be chosen from the 2 Ethicists and Cassie.
    
    # The 5 people not acting as direct neighbors can be arranged in 5! ways in the remaining seats.
    ways_remaining_people = math.factorial(5)

    # Case A: The neighbors of the Super-Block are the 2 Ethicists.
    # 2! ways to arrange the ethicists. Internal arrangements can be any of the N_total.
    ways_ethicist_neighbors = 2 * ways_remaining_people * N_total
    
    # Case B: The neighbors are Cassie and one of the Ethicists.
    # Cassie can be at the S-end or M-end.
    # There are 2 choices for the Ethicist.
    # If Cassie is at S-end, internal arrangements must be N_S_female.
    # If Cassie is at M-end, internal arrangements must be N_M_female.
    ways_cassie_neighbor = (2 * ways_remaining_people * N_S_female) + \
                           (2 * ways_remaining_people * N_M_female)

    # Total number of arrangements is the sum of all valid cases.
    total_ways = ways_ethicist_neighbors + ways_cassie_neighbor

    # Print the step-by-step calculation
    print("Step 1: Calculate the number of internal arrangements for the Scientist-Mathematician super-block.")
    print(f"N_FF (Female Scientist end, Female Mathematician end): {N_FF}")
    print(f"N_FM (Female Scientist end, Male Mathematician end):   {N_FM}")
    print(f"N_MF (Male Scientist end, Female Mathematician end):   {N_MF}")
    print(f"N_MM (Male Scientist end, Male Mathematician end):     {N_MM}")
    print("-" * 30)
    print(f"Total internal arrangements (N_total): {N_total}")
    print(f"Internal arrangements with Female Scientist end (N_S_female): {N_S_female}")
    print(f"Internal arrangements with Female Mathematician end (N_M_female): {N_M_female}")
    print("\nStep 2: Calculate arrangements based on the neighbors of the super-block.")
    print(f"Ways to arrange the 5 non-neighboring people: 5! = {ways_remaining_people}")
    print(f"Ways if neighbors are Ethicists = 2 * {ways_remaining_people} * {N_total} = {ways_ethicist_neighbors}")
    print(f"Ways if Cassie is a neighbor = (2 * {ways_remaining_people} * {N_S_female}) + (2 * {ways_remaining_people} * {N_M_female}) = {ways_cassie_neighbor}")
    print("\nStep 3: Calculate the final total.")
    print(f"Total Ways = {ways_ethicist_neighbors} + {ways_cassie_neighbor}")
    print(f"Total Ways = {total_ways}")

    return total_ways

# Run the calculation and print the final answer.
final_answer = solve_seating_arrangement()
print(f"\nFinal Answer: The total number of ways to arrange the table is {final_answer}.")
print(f"<<<{final_answer}>>>")
