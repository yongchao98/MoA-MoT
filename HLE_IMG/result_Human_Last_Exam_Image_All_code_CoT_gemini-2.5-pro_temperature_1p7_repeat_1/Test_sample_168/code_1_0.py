import math

def solve_light_path_problem():
    """
    Calculates the number of possible light ray paths based on reflection constraints.
    """

    # --- Introduction ---
    print("This program calculates the number of ways a light ray can travel from point M to N.")
    print("The ray must reflect twice on mirror G1 and once on mirrors G2, G3, and G4.")
    print("The multiset of reflections is {G1, G1, G2, G3, G4}.\n")

    print("A key physical constraint is that a ray cannot reflect on the same mirror twice consecutively.")
    print("This means the two G1 reflections cannot be adjacent in the sequence.\n")

    print("Our method is to find the total permutations and subtract the forbidden ones where G1s are together.\n")

    # --- Step 1: Calculate Total Permutations ---
    total_mirrors = 5
    repeated_mirrors = 2
    
    print(f"Step 1: Calculate the total number of permutations for the {total_mirrors} reflections.")
    print("This is the permutation of the multiset {G1, G1, G2, G3, G4}.")
    
    # Formula for permutations with repetition: n! / k!
    total_permutations = math.factorial(total_mirrors) / math.factorial(repeated_mirrors)
    
    print(f"Equation for total permutations: {total_mirrors}! / {repeated_mirrors}!")
    print(f"Calculation: {math.factorial(total_mirrors)} / {math.factorial(repeated_mirrors)} = {int(total_permutations)}\n")

    # --- Step 2: Calculate Forbidden Permutations ---
    # We treat (G1, G1) as a single unit. Now we permute {(G1,G1), G2, G3, G4}.
    items_if_adjacent = 4
    
    print(f"Step 2: Calculate the number of 'forbidden' permutations where the two G1s are adjacent.")
    print("We can treat the pair (G1,G1) as a single, unique item.")
    print(f"This reduces the problem to finding the permutations of {items_if_adjacent} distinct items: {(G1,G1), G2, G3, G4}.")
    
    forbidden_permutations = math.factorial(items_if_adjacent)
    
    print(f"Equation for forbidden permutations: {items_if_adjacent}!")
    print(f"Calculation: {forbidden_permutations}\n")

    # --- Step 3: Final Result ---
    number_of_ways = total_permutations - forbidden_permutations
    
    print("Step 3: Subtract the forbidden permutations from the total to find the number of valid paths.")
    print(f"Final Equation: Total Permutations - Forbidden Permutations")
    print(f"Result: {int(total_permutations)} - {forbidden_permutations} = {int(number_of_ways)}")

solve_light_path_problem()
<<<36>>>