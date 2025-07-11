import math

# The problem is to find the number of valid sequences for a light ray reflecting
# twice on mirror G1, and once each on mirrors G2, G3, and G4.

# Step 1: Calculate the total number of theoretical reflection sequences.
# The set of reflections is {G1, G1, G2, G3, G4}.
# This is a permutation with repetition problem for 5 items, where one item (G1) is repeated twice.
# The total number of permutations is given by the formula n! / k!, where n is the total number
# of reflections and k is the number of times an item is repeated.
total_reflections = 5
repeated_reflections = 2
total_sequences = math.factorial(total_reflections) / math.factorial(repeated_reflections)

# Step 2: Calculate the number of impossible sequences.
# A sequence is physically impossible if the ray reflects off the same mirror twice in a row.
# We need to count the sequences where the two G1 reflections are adjacent.
# We can do this by treating the '(G1, G1)' pair as a single block.
# Now, we are permuting 4 distinct items: {(G1,G1), G2, G3, G4}.
# The number of permutations of these 4 items is 4!.
items_for_invalid_sequence = 4
impossible_sequences = math.factorial(items_for_invalid_sequence)

# Step 3: Calculate the number of valid ways.
# This is the total number of sequences minus the number of impossible sequences.
valid_ways = total_sequences - impossible_sequences

# Print the step-by-step calculation.
print("To find the number of ways, we follow these steps:")
print("1. Calculate the total number of ways to arrange the 5 reflections {G1, G1, G2, G3, G4}:")
print(f"   Total arrangements = 5! / 2! = {int(math.factorial(5))} / {int(math.factorial(2))} = {int(total_sequences)}")
print("\n2. Calculate the number of impossible arrangements where the two G1 reflections are consecutive.")
print("   Treat '(G1, G1)' as a single block and arrange {(G1,G1), G2, G3, G4}:")
print(f"   Impossible arrangements = 4! = {int(impossible_sequences)}")
print("\n3. Subtract the impossible arrangements from the total to get the number of valid ways:")
print("   Valid ways = Total arrangements - Impossible arrangements")
print(f"   Final Equation: {int(total_sequences)} - {int(impossible_sequences)} = {int(valid_ways)}")
