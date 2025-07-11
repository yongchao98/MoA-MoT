import math

# Plan Explanation
print("This problem asks for the number of ways a light ray path can be formed.")
print("The path involves 2 reflections on mirror G1, and 1 reflection each on G2, G3, and G4.")
print("A physically valid path cannot have consecutive reflections on the same mirror.")
print("Therefore, we need to count the number of sequences of reflections where the two G1s are not adjacent.")
print("-" * 20)

# Step 1: Arrange the unique mirrors (G2, G3, G4).
# There are 3 distinct mirrors that are each hit once.
num_unique_mirrors = 3
perms_unique_mirrors = math.factorial(num_unique_mirrors)
print(f"First, we calculate the number of ways to arrange the mirrors G2, G3, and G4.")
print(f"Number of arrangements for 3 mirrors = 3! = {perms_unique_mirrors}")
print()

# Step 2: Create slots for the G1 reflections.
# Arranging 3 mirrors creates 4 slots where other items can be placed so they are not adjacent.
# Example: _ G2 _ G3 _ G4 _
num_slots = num_unique_mirrors + 1
print(f"An arrangement of {num_unique_mirrors} mirrors creates {num_slots} slots to place the G1 reflections.")

# Step 3: Place the two G1 reflections into the available slots.
# We need to choose 2 distinct slots from the 4 available ones.
num_g1_reflections = 2
ways_to_place_g1s = math.comb(num_slots, num_g1_reflections)
print(f"We need to place {num_g1_reflections} 'G1' reflections into these {num_slots} slots.")
print(f"The number of ways to choose 2 slots from 4 is C(4, 2).")
print(f"C(4, 2) = 4! / (2! * (4-2)!) = {ways_to_place_g1s}")
print()

# Step 4: Calculate the total number of valid sequences.
total_ways = perms_unique_mirrors * ways_to_place_g1s
print("The total number of ways is the product of these two results.")
print("Final Equation:")
print(f"Total ways = (permutations of unique mirrors) * (ways to place G1s)")
print(f"Total ways = {perms_unique_mirrors} * {ways_to_place_g1s} = {total_ways}")