import math

# Step 1: Define the geometric parameters based on the analysis.
# Radius of a large white circle (R)
R = 5
# Radius of a small yellow circle (r_y)
r_y = 1
# Number of yellow circles
num_yellow_circles = 6

# Step 2: Calculate the total dimensions of the target.
# Width of the green rectangle (w_g), derived from integer constraints
w_g = 8
# Height of the target (H), from H = 4*R + 2
H = 4 * R + 2
# Width of the target (W), derived from the layout of the middle row
# W = w_g (left rect) + R (to 1st center) + 2R (to 2nd center) + 2R (to 3rd center) + R+r_y (to yellow center) + r_y (to edge)
W = w_g + 5 * R + r_y + (R + r_y) # A more detailed calculation based on centers
W = 8 + 5*5 + 1 + (5+1) # Simplified: W = 40 from analysis
W = 40

# Step 3: Calculate the total areas.
# Total area of the target
total_area = W * H
# Total area of the yellow circles
area_one_yellow = math.pi * r_y**2
total_yellow_area = num_yellow_circles * area_one_yellow

# Step 4: Calculate the expected number of hits for 10000 shots.
num_shots = 10000
expected_hits = (total_yellow_area / total_area) * num_shots

# Step 5: Print the explanation and the final equation with all numbers.
print("To find the expected number of hits, we calculate the ratio of the yellow area to the total area and multiply by the number of shots.")
print("\n--- Geometric Parameters ---")
print(f"Radius of a large white circle (R): {R} cm")
print(f"Radius of a small yellow circle (r_y): {r_y} cm")
print(f"Number of yellow circles: {num_yellow_circles}")

print("\n--- Target Dimensions ---")
print(f"Target Width (W): {W} cm")
print(f"Target Height (H): {H} cm")

print("\n--- Area Calculations ---")
print(f"Area of one yellow circle = pi * ({r_y})^2 = {area_one_yellow:.4f} cm^2")
print(f"Total area of yellow circles = {num_yellow_circles} * {area_one_yellow:.4f} = {total_yellow_area:.4f} cm^2")
print(f"Total area of the target = {W} * {H} = {total_area} cm^2")

print("\n--- Final Calculation ---")
print("Expected Hits = (Total Yellow Area / Total Target Area) * Number of Shots")
print(f"Expected Hits = (({num_yellow_circles} * pi * {r_y}**2) / ({W} * {H})) * {num_shots}")
print(f"Expected Hits = {expected_hits:.4f}")

# Final Answer
final_answer = (6 * math.pi / 880) * 10000
print(f"\nFinal numerical answer: {final_answer}")
<<<214.1995>>>