import math

def solve_mirror_problem():
    """
    Calculates the number of ways a light ray can reflect according to the problem's rules.
    """
    # The problem is to find the number of sequences of reflections {G1, G1, G2, G3, G4}
    # such that the two G1 reflections are not consecutive.

    # Step 1: Arrange the mirrors other than G1. These are G2, G3, and G4.
    num_other_mirrors = 3
    # The number of ways to arrange these 3 distinct mirrors is 3!
    perms_other_mirrors = math.factorial(num_other_mirrors)

    # Step 2: These arrangements create slots where the G1 mirrors can be placed.
    # For an arrangement of 3 mirrors, there are 4 slots (_M1_M2_M3_).
    num_slots = num_other_mirrors + 1

    # Step 3: We need to place the 2 G1 mirrors into these slots.
    num_g1_mirrors = 2
    # The number of ways to choose 2 slots from 4 is given by the combination formula C(n, k).
    ways_to_place_g1s = math.comb(num_slots, num_g1_mirrors)

    # Step 4: The total number of ways is the product of the permutations and combinations.
    total_ways = perms_other_mirrors * ways_to_place_g1s
    
    # Print the explanation and the final equation.
    print(f"The number of ways is calculated by arranging the 3 distinct mirrors (G2, G3, G4) and then placing the 2 identical mirrors (G1) in the available slots.")
    print(f"Number of ways to arrange G2, G3, G4 = {num_other_mirrors}! = {perms_other_mirrors}")
    print(f"Number of slots created = {num_slots}")
    print(f"Number of ways to place the {num_g1_mirrors} G1 mirrors in {num_slots} slots = C({num_slots}, {num_g1_mirrors}) = {ways_to_place_g1s}")
    print("\nFinal Equation:")
    print(f"Total Number of Ways = {num_other_mirrors}! * C({num_slots}, {num_g1_mirrors}) = {perms_other_mirrors} * {ways_to_place_g1s} = {total_ways}")

solve_mirror_problem()
<<<36>>>