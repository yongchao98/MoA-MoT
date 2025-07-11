import math

# Based on the geometric analysis of the problem's constraints:
# - The radius of a yellow circle is given as r = 1 cm.
# - The radius of a white circle is derived to be R = 4 cm.
# - The total width of the target is derived to be W = 32 cm.
# - The total height of the target is derived to be H = 20 cm.

# Parameters
r_yellow = 1
num_yellow_circles = 6
total_shots = 10000
W = 32
H = 20

# Calculate the total area of the target
total_area = W * H

# Calculate the total area of the yellow circles
area_one_yellow_circle = math.pi * r_yellow**2
total_yellow_area = num_yellow_circles * area_one_yellow_circle

# Calculate the expected number of hits
expected_hits = total_shots * (total_yellow_area / total_area)

# Print the final equation with all the numbers
print("Expected Hits = Total Shots * (Total Yellow Area / Total Target Area)")
print(f"Expected Hits = {total_shots} * ({num_yellow_circles} * pi * {r_yellow}^2) / ({W} * {H})")
print(f"Expected Hits = {total_shots} * ({num_yellow_circles:.1f} * {math.pi:.4f}) / {total_area}")
print(f"Expected Hits = {expected_hits:.2f}")

# The number of hits must be an integer, so we round the result.
final_answer = round(expected_hits)
print(f"\nRounding to the nearest integer, the expected number of hits is {final_answer}.")
<<<295>>>