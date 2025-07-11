import math

# Step 1: Define the number of reflections on each type of mirror.
# Vertical mirrors are G1 and G3. Reflections: 2 on G1, 1 on G3. Total = 3.
num_v_reflections = 3
# The multiset of vertical mirrors is {G1, G1, G3}.
# We have 2 identical G1 mirrors.
count_g1 = 2
count_g3 = 1

# Horizontal mirrors are G2 and G4. Reflections: 1 on G2, 1 on G4. Total = 2.
num_h_reflections = 2
# The set of horizontal mirrors is {G2, G4}. All distinct.
count_g2 = 1
count_g4 = 1

# Step 2: Calculate the number of ways to arrange the vertical mirrors.
# The formula for permutations of a multiset is n! / (n1! * n2! * ...).
# For vertical mirrors {G1, G1, G3}: 3! / 2!
num_v_arrangements = math.factorial(num_v_reflections) // math.factorial(count_g1)
print(f"The number of reflections on vertical mirrors {G1, G3} is 3.")
print(f"The number of ways to arrange the vertical mirrors is {num_v_reflections}! / {count_g1}! = {num_v_arrangements}")
print("-" * 20)

# Step 3: Calculate the number of ways to arrange the horizontal mirrors.
# For horizontal mirrors {G2, G4}: 2! / (1! * 1!) which is just 2!
num_h_arrangements = math.factorial(num_h_reflections)
print(f"The number of reflections on horizontal mirrors {G2, G4} is 2.")
print(f"The number of ways to arrange the horizontal mirrors is {num_h_reflections}! = {num_h_arrangements}")
print("-" * 20)

# Step 4: Calculate the total number of ways.
# The total number of unique paths is the product of the two arrangements.
total_ways = num_v_arrangements * num_h_arrangements
print("The physical constraints require the reflection sequence to alternate between vertical and horizontal mirrors.")
print("With 3 vertical and 2 horizontal reflections, the pattern must be V-H-V-H-V.")
print("The total number of ways to draw the light ray path is the product of the possible arrangements for each type.")
print(f"Total ways = (Ways for Vertical) * (Ways for Horizontal)")
print(f"Total ways = {num_v_arrangements} * {num_h_arrangements} = {total_ways}")

# Final Answer format
print("\n<<<6>>>")
