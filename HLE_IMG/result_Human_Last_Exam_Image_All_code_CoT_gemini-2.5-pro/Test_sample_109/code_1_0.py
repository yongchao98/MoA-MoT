import math

# Step 1: Define the geometric parameters based on the analysis
# Radius of yellow circles
r_yellow = 1
# Radius of white circles, derived from geometric constraints
R_white = 5
# Height of the target
H = 26
# Width of the target
W = 32

# Step 2: Calculate the total area of the yellow circles
num_yellow_circles = 6
area_per_yellow_circle = math.pi * r_yellow**2
total_yellow_area = num_yellow_circles * area_per_yellow_circle

# Step 3: Calculate the total area of the target
total_target_area = W * H

# Step 4: Calculate the probability of hitting a yellow circle
probability_hit_yellow = total_yellow_area / total_target_area

# Step 5: Calculate the expected number of hits in 10000 shots
total_shots = 10000
expected_hits = total_shots * probability_hit_yellow

# Print the final equation with all the numbers
print(f"Expected Hits = Total Shots * (Total Yellow Area / Total Target Area)")
print(f"Expected Hits = {total_shots} * ({num_yellow_circles} * pi * {r_yellow}^2) / ({W} * {H})")
print(f"Expected Hits = {total_shots} * ({num_yellow_circles * area_per_yellow_circle:.4f}) / ({total_target_area})")
print(f"Expected Hits = {expected_hits:.4f}")

# The question asks for "how many times", which implies an integer or a decimal value for an expectation.
# We will provide the decimal value, as is standard for expected values.
final_answer = round(expected_hits, 4)