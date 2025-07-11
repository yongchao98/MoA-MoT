import math

def solve_light_path_ways():
    """
    This function calculates the number of ways a light ray can travel from M to N
    with 2 reflections on G1, and 1 on G2, G3, and G4.
    """
    
    # Step 1: Define the parameters of the multiset of reflections {G1, G1, G2, G3, G4}.
    total_reflections = 5
    g1_reflections = 2
    
    print("Step 1: Calculate the total number of permutations of the 5 reflections {G1, G1, G2, G3, G4}.")
    # The formula for permutations of a multiset is n! / (n1! * n2! * ...).
    # Here, n=5 total reflections, with n1=2 identical reflections on G1.
    total_permutations = math.factorial(total_reflections) // math.factorial(g1_reflections)
    print(f"Total permutations = {total_reflections}! / {g1_reflections}! = {math.factorial(total_reflections)} / {math.factorial(g1_reflections)} = {total_permutations}")
    print("-" * 40)
    
    # Step 2: Calculate the number of forbidden permutations where the two G1 reflections are adjacent.
    # To do this, we treat the two G1s as a single, combined item: '(G1,G1)'.
    # Now, we are permuting a set of 4 distinct items: {(G1,G1), G2, G3, G4}.
    print("Step 2: Calculate the number of forbidden paths where the two G1 reflections are consecutive.")
    print("We treat '(G1,G1)' as a single block, and permute it with G2, G3, G4.")
    num_items_for_forbidden_case = (total_reflections - g1_reflections) + 1
    forbidden_permutations = math.factorial(num_items_for_forbidden_case)
    print(f"Number of forbidden permutations = {num_items_for_forbidden_case}! = {forbidden_permutations}")
    print("-" * 40)

    # Step 3: Subtract the forbidden permutations from the total to get the number of valid ways.
    valid_ways = total_permutations - forbidden_permutations
    print("Step 3: The final number of ways is the total permutations minus the forbidden permutations.")
    print(f"Final Calculation: {total_permutations} - {forbidden_permutations} = {valid_ways}")
    print("-" * 40)
    
    print("The final equation with all numbers is:")
    print(f"({math.factorial(total_reflections)} / {math.factorial(g1_reflections)}) - {math.factorial(num_items_for_forbidden_case)} = {valid_ways}")

solve_light_path_ways()
<<<36>>>