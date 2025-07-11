import math

# Step 1 & 2: Define the determined dimensions based on geometric analysis.
r_yellow = 1
R_white = 40  # Determined by finding a valid Pythagorean triple solution.

# Step 3: Calculate the total area of the target.
W_target = 8 * R_white
H_target = 148  # Determined from the Pythagorean triple and other constraints.
total_area = W_target * H_target

# Step 4: Calculate the total area of the yellow circles.
num_yellow_circles = 6  # By counting from the image.
area_one_yellow = math.pi * (r_yellow ** 2)
total_yellow_area = num_yellow_circles * area_one_yellow

# Step 5: Calculate the expected number of hits in 10000 shots.
total_shots = 10000
expected_hits = total_shots * (total_yellow_area / total_area)

print(f"The radius of a large white circle (R) is {R_white} cm.")
print(f"The radius of a yellow circle (r) is {r_yellow} cm.")
print(f"The width of the target is 8 * R = 8 * {R_white} = {W_target} cm.")
print(f"The height of the target is {H_target} cm.")
print(f"The total area of the target is {W_target} * {H_target} = {total_area} cm^2.")
print(f"The number of yellow circles is {num_yellow_circles}.")
print(f"The area of one yellow circle is pi * {r_yellow}^2 = {area_one_yellow:.4f} cm^2.")
print(f"The total area of yellow circles is {num_yellow_circles} * pi = {total_yellow_area:.4f} cm^2.")
print("\nThe expected number of hits in 10000 shots is calculated as:")
print(f"Expected Hits = {total_shots} * (Total Yellow Area / Total Target Area)")
print(f"Expected Hits = {total_shots} * ({num_yellow_circles} * pi / {total_area})")
print(f"Expected Hits = {expected_hits:.4f}")
# The question asks for "how many times", which suggests a numerical value.
# Let's provide the rounded integer value as well as the more precise float.
print(f"\nRounded to the nearest whole number, you can expect to hit a yellow circle about {round(expected_hits)} times.")
