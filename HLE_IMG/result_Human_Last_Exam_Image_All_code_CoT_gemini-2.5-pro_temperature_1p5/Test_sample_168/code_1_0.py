import math

def solve_light_path_ways():
    """
    Calculates the number of ways a light ray can travel reflecting
    twice on G1 and once on G2, G3, and G4.
    """
    print("This script calculates the number of possible light ray paths based on combinatorial principles.")
    print("------------------------------------------------------------------------------------------\n")

    # Step 1: Define the multiset of reflections.
    # The unique mirrors are G2, G3, G4. The repeating mirror is G1 (twice).
    num_unique_mirrors = 3
    num_repeating_mirrors = 2

    # Step 2: Arrange the unique mirrors {G2, G3, G4}.
    # The number of permutations is 3!
    perms_of_unique_mirrors = math.factorial(num_unique_mirrors)
    print(f"Step 1: Calculate the number of ways to arrange the unique mirrors (G2, G3, G4).")
    print(f"This is 3! = {perms_of_unique_mirrors}.\n")

    # Step 3: For each arrangement, determine the valid ways to insert the two G1 mirrors.
    # An arrangement of 3 mirrors creates 4 potential slots for insertion (e.g., _ G2 _ G3 _ G4 _).
    num_slots = num_unique_mirrors + 1
    print(f"Step 2: For each arrangement, find the number of ways to insert the two G1 mirrors.")
    print(f"An arrangement like 'G2-G3-G4' creates {num_slots} slots: _ G2 _ G3 _ G4 _\n")

    # To avoid the 'G1 G1' pattern, the two G1s must be in different slots.
    # The number of ways to choose 2 distinct slots from 4 is C(4, 2).
    total_placements = math.comb(num_slots, num_repeating_mirrors)
    print(f"  a) To avoid the 'G1-G1' pattern, we place the two G1s in {num_repeating_mirrors} different slots.")
    print(f"     Number of ways to choose {num_repeating_mirrors} slots from {num_slots} = C({num_slots}, {num_repeating_mirrors}) = {total_placements}.\n")

    # From these placements, we must exclude those forming forbidden 'G1-X-G1' patterns.
    # The mirrors perpendicular to G1 are G2 and G4. So, 'G1-G2-G1' and 'G1-G4-G1' are forbidden.
    # For any permutation of {G2, G3, G4}, there is exactly one way to form 'G1-G2-G1' (by placing G1s around G2)
    # and exactly one way to form 'G1-G4-G1' (by placing G1s around G4).
    num_forbidden_adj_patterns = 2
    print(f"  b) We must exclude placements creating 'G1-G2-G1' or 'G1-G4-G1' patterns.")
    print(f"     For each arrangement, there are {num_forbidden_adj_patterns} such forbidden placements.\n")
    
    # Calculate the number of valid placements per arrangement.
    valid_placements_per_perm = total_placements - num_forbidden_adj_patterns
    print(f"  c) The number of valid placements for each arrangement is {total_placements} - {num_forbidden_adj_patterns} = {valid_placements_per_perm}.\n")

    # Step 4: Calculate the total number of ways.
    total_ways = perms_of_unique_mirrors * valid_placements_per_perm
    print("Step 3: Calculate the total number of ways.")
    print("Total ways = (Arrangements of G2,G3,G4) * (Valid placements of G1s per arrangement)\n")
    
    # The final equation and result
    print("Final Equation:")
    print(f"Number of ways = {perms_of_unique_mirrors} * ({total_placements} - {num_forbidden_adj_patterns})")
    print(f"Number of ways = {perms_of_unique_mirrors} * {valid_placements_per_perm}")
    print(f"Number of ways = {total_ways}")

solve_light_path_ways()
<<<24>>>