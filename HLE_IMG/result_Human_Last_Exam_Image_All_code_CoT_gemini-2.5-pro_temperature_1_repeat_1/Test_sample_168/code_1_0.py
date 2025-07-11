import math

def solve_light_path_ways():
    """
    Calculates the number of ways a light ray can reflect off a set of mirrors
    based on combinatorial principles.
    """

    # Step 1: Define the multiset of reflections.
    # Total reflections = 5, with G1 appearing twice.
    total_reflections = 5
    g1_reflections = 2
    other_mirrors = ["G2", "G3", "G4"]
    num_other_mirrors = len(other_mirrors)

    print("Step 1: The multiset of reflections is {G1, G1, G2, G3, G4}.")
    print(f"There are a total of {total_reflections} reflections.")
    print("-" * 40)

    # Step 2: Calculate the total number of permutations of this multiset.
    # The formula for permutations of a multiset is N! / (n1! * n2! * ...).
    # Here, N=5 and n1=2 for G1.
    total_sequences = math.factorial(total_reflections) / math.factorial(g1_reflections)
    
    print("Step 2: Calculate the total number of possible sequences (permutations).")
    print(f"Total sequences = {total_reflections}! / {g1_reflections}! = {int(total_sequences)}")
    print("-" * 40)

    # Step 3: Identify and calculate the number of invalid sequences.
    # An invalid sequence is one where the two G1 reflections are consecutive.
    # We treat '(G1, G1)' as a single block.
    # We now have to permute 4 items: {(G1, G1), G2, G3, G4}.
    items_in_invalid_case = num_other_mirrors + 1  # The block (G1,G1) + other mirrors
    invalid_sequences = math.factorial(items_in_invalid_case)

    print("Step 3: Calculate the number of invalid sequences (where G1s are adjacent).")
    print("We treat the two G1 reflections as a single block.")
    print(f"The number of items to permute becomes {items_in_invalid_case}.")
    print(f"Number of invalid sequences = {items_in_invalid_case}! = {invalid_sequences}")
    print("-" * 40)

    # Step 4: Calculate the final number of valid ways.
    # Valid ways = Total sequences - Invalid sequences.
    valid_ways = total_sequences - invalid_sequences
    
    print("Step 4: Calculate the number of valid ways.")
    print("Valid ways = Total sequences - Invalid sequences")
    print(f"The final number of ways is {int(total_sequences)} - {invalid_sequences} = {int(valid_ways)}")

solve_light_path_ways()
