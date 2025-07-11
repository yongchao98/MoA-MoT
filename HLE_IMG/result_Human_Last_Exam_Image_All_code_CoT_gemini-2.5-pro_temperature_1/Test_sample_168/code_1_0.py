import math

def calculate_reflection_ways():
    """
    Calculates the number of ways a light ray can reflect off a set of mirrors
    based on combinatorial principles.
    """
    # Total number of reflections in the sequence.
    n_total = 5  # {G1, G1, G2, G3, G4}

    # Number of reflections on the repeated mirror (G1).
    n_repeated = 2

    # Step 1: Calculate the total number of unique permutations of the multiset.
    # Formula: n! / k!, where n is total items, and k is the count of repeated items.
    # Here, it is 5! / 2!
    total_permutations = math.factorial(n_total) // math.factorial(n_repeated)

    # Step 2: Calculate the number of invalid permutations where the two G1s are adjacent.
    # We treat the two G1s as a single block, so we are permuting 4 distinct items:
    # {(G1,G1), G2, G3, G4}.
    # The number of items to permute becomes (n_total - n_repeated + 1) = 4.
    n_items_for_invalid_perms = n_total - n_repeated + 1
    invalid_permutations = math.factorial(n_items_for_invalid_perms)

    # Step 3: The number of valid ways is the total minus the invalid permutations.
    valid_ways = total_permutations - invalid_permutations

    # Step 4: Print the final equation with the calculated numbers.
    print(f"{total_permutations} - {invalid_permutations} = {valid_ways}")

calculate_reflection_ways()