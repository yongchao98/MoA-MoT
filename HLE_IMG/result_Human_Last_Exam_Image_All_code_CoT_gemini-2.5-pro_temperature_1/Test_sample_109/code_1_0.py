import math

# --- Step 1: Define the dimensions based on the problem description ---
# Radius of a yellow circle (given)
r_y = 1
# Radius of a large white circle, derived from R = 4 * r_y
R = 4

# --- Step 2: Determine the target dimensions ---
# Height of the target, derived from integer coordinate constraints
H = 24
# Width of the target, determined by the widest row of circles (3 large circles)
W = 3 * (2 * R)

# --- Step 3: Calculate Areas ---
# Total area of the target
total_area = W * H
# Number of yellow circles on the target
num_yellow_circles = 6
# Area of a single yellow circle
yellow_circle_area = math.pi * r_y**2
# Total area of all yellow circles
total_yellow_area = num_yellow_circles * yellow_circle_area

# --- Step 4: Calculate the expected number of hits ---
# Number of shots
num_shots = 10000
# Probability of hitting a yellow circle
probability_hit_yellow = total_yellow_area / total_area
# Expected number of hits
expected_hits = num_shots * probability_hit_yellow

# --- Print the logic and the final result ---
print("--- Calculation Steps ---")
print(f"Radius of a yellow circle (r_y): {r_y} cm")
print(f"Radius of a large white circle (R): {R} cm")
print(f"Total height of the target (H): {H} cm")
print(f"Total width of the target (W): {W} cm")
print(f"Total area of the target: {W} * {H} = {total_area} cm^2")
print(f"Number of yellow circles: {num_yellow_circles}")
print(f"Area of one yellow circle: pi * {r_y}^2 = {yellow_circle_area:.4f} cm^2")
print(f"Total area of yellow circles: {num_yellow_circles} * pi = {total_yellow_area:.4f} cm^2")
print(f"Probability of hitting a yellow circle: ({num_yellow_circles} * pi) / {total_area} = {probability_hit_yellow:.6f}")
print("\n--- Final Calculation ---")
print(f"Expected number of hits in {num_shots} shots is:")
print(f"{num_shots} * ({num_yellow_circles} * pi / {total_area})")
print(f"= {num_shots} * {probability_hit_yellow:.6f}")
print(f"= {expected_hits:.4f}")
print(f"\nThe final answer is {expected_hits}")
