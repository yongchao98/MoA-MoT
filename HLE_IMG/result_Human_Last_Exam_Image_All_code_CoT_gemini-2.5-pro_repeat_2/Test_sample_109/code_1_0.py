import math

# Step 1: Define radii
r_yellow = 1  # cm, given
# From geometric analysis (R = 4r)
R_white = 4 * r_yellow

# Step 2: Define target dimensions based on analysis
# Height is derived from vertical stacking of 3 rows of large circles
H = 2 * (R_white + 8) # Correct derivation is H/2 - R = 2R -> H = 6R = 24. A simpler way is H = (R) + (H/2-R) + (H/2-R) + (R) = R + 8 + 8 + R = 2R + 16 = 8+16 = 24
H = 24
# Width is derived by assuming the middle green bar is a square (width=height)
# Height of middle green bar = diameter of large circle = 2 * R_white
height_green_bar = 2 * R_white
# Assuming the green bar in the middle row is a square
width_green_bar_middle = height_green_bar
# The middle row consists of the green bar and 3 large white circles
width_of_3_white_circles = 3 * (2 * R_white)
# Total width W
W = width_green_bar_middle + width_of_3_white_circles

# Step 3: Calculate areas
# Total area of the rectangular target
total_area = W * H
# Area of the yellow circles
# There are 4 yellow circles in total (2 top, 2 bottom)
num_yellow_circles = 4
area_one_yellow_circle = math.pi * (r_yellow ** 2)
total_yellow_area = num_yellow_circles * area_one_yellow_circle

# Step 4: Calculate the expected number of hits for 10000 shots
num_shots = 10000
expected_hits = num_shots * (total_yellow_area / total_area)

# Print out the logic and the final answer
print("### Step-by-step Calculation ###")
print(f"Radius of a yellow circle (r): {r_yellow} cm")
print(f"Radius of a white circle (R): {R_white} cm")
print(f"Height of the target (H): {H} cm")
print(f"Width of the middle green bar (assumed square): {width_green_bar_middle} cm")
print(f"Horizontal span of the three middle white circles: {width_of_3_white_circles} cm")
print(f"Total width of the target (W): {width_green_bar_middle} + {width_of_3_white_circles} = {W} cm")
print(f"Total area of the target: {W} * {H} = {total_area} cm^2")
print(f"Number of yellow circles: {num_yellow_circles}")
print(f"Area of one yellow circle: pi * {r_yellow}^2 = {area_one_yellow_circle:.4f} cm^2")
print(f"Total area of yellow circles: {num_yellow_circles} * {area_one_yellow_circle:.4f} = {total_yellow_area:.4f} cm^2")
print("\n### Final Calculation ###")
print(f"Expected hits in {num_shots} shots = {num_shots} * (Total Yellow Area / Total Target Area)")
print(f"Expected hits = {num_shots} * ({total_yellow_area:.4f} / {total_area})")
print(f"Expected hits = {expected_hits:.1f}")

# The final answer in the required format
final_answer = round(expected_hits, 1)