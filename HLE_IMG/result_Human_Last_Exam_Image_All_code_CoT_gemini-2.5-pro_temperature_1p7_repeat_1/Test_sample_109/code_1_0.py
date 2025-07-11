import math

# Step 1 & 2: Define constants based on the geometric analysis.
r = 1  # Radius of yellow circles in cm
R = 4  # Radius of white circles in cm

# Step 3, 4, 5: Determine the height of the target.
H = 28  # Total height of the target in cm

# Step 6, 7, 8: Determine the width of the target.
w_g = 8  # Assumed width of the green rectangle, equal to the diameter of a white circle.
num_white_circles_in_row = 4
W = w_g + num_white_circles_in_row * (2 * R)

# Step 9: Determine the number of yellow circles.
# Based on the model, there are 3 valleys in the top row and 3 in the bottom row containing yellow circles.
num_yellow_circles = 6

# Step 10 & 11: Calculate the total area of yellow circles and the total area of the target.
area_one_yellow_circle = math.pi * r**2
total_yellow_area = num_yellow_circles * area_one_yellow_circle
total_target_area = W * H

# Calculate the probability of hitting a yellow circle.
probability = total_yellow_area / total_target_area

# Calculate the expected number of hits in 10000 shots.
num_shots = 10000
expected_hits = num_shots * probability

print("Step-by-step calculation:")
print(f"Radius of a small yellow circle (r): {r} cm")
print(f"Radius of a large white circle (R): {R} cm")
print(f"Total height of the target (H): {H} cm")
print(f"Total width of the target (W): {W} cm")
print(f"Number of yellow circles (N_y): {num_yellow_circles}")
print(f"Area of one yellow circle: pi * {r}^2 = {area_one_yellow_circle:.4f} cm^2")
print(f"Total area of yellow circles: {num_yellow_circles} * pi * {r}^2 = {total_yellow_area:.4f} cm^2")
print(f"Total area of the target: {W} * {H} = {total_target_area} cm^2")
print("\nFinal Equation for Expected Hits:")
# To satisfy the output format instruction "output each number in the final equation"
# E = num_shots * (num_yellow_circles * pi * r^2) / (W * H)
print(f"Expected Hits = {num_shots} * ({num_yellow_circles} * pi * {r}^2) / ({W} * {H})")
print(f"Expected Hits = {num_shots} * {num_yellow_circles} / ({W} * {H}) * pi")
print(f"Expected Hits = {expected_hits:.1f}")
