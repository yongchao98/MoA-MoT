import math

def solve_mirror_problem():
    """
    Calculates the number of ways a light ray can travel reflecting
    twice on mirror G1 and once on mirrors G2, G3, and G4.
    """
    
    # Step 1: Define the number of reflections.
    m = 2  # Reflections on G1
    num_other_mirrors = 3 # G2, G3, G4
    total_reflections = m + num_other_mirrors
    
    print("This problem is equivalent to finding the number of valid sequences of reflections.")
    print(f"The set of reflections is {{G1, G1, G2, G3, G4}}.")
    print(f"Total number of reflections = {m} (on G1) + 1 (on G2) + 1 (on G3) + 1 (on G4) = {total_reflections}\n")

    # Step 2: Calculate the total number of permutations of the multiset.
    # Formula: n! / m! where n is total reflections and m is the count of the repeated reflection.
    total_permutations_numerator = math.factorial(total_reflections)
    total_permutations_denominator = math.factorial(m)
    total_permutations = total_permutations_numerator // total_permutations_denominator

    print("Step 1: Calculate the total number of permutations for the multiset of reflections.")
    print(f"Total permutations = {total_reflections}! / {m}! = {total_permutations_numerator} / {total_permutations_denominator} = {total_permutations}\n")

    # Step 3: Identify and calculate the number of invalid permutations.
    # An invalid permutation is one where the two G1 reflections are adjacent (e.g., ...G1, G1...).
    # To count these, we treat the two G1s as a single block.
    # We are now permuting (n - m + 1) items.
    items_for_invalid_perm = (total_reflections - m + 1)
    invalid_permutations = math.factorial(items_for_invalid_perm)

    print("Step 2: A light ray cannot reflect on the same mirror twice in a row.")
    print("We must subtract the number of invalid sequences where the two G1 reflections are adjacent.")
    print("To count these, we treat '(G1, G1)' as a single block.")
    print(f"The number of items to permute becomes {items_for_invalid_perm} (i.e., {{(G1,G1), G2, G3, G4}}).")
    print(f"Number of invalid permutations = {items_for_invalid_perm}! = {invalid_permutations}\n")

    # Step 4: Calculate the final number of valid ways.
    valid_ways = total_permutations - invalid_permutations

    print("Step 3: The number of valid ways is the total permutations minus the invalid permutations.")
    print(f"Number of ways = {total_permutations} - {invalid_permutations} = {valid_ways}")

    # Final answer in the required format
    print(f"\n<<<{valid_ways}>>>")

solve_mirror_problem()