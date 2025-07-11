import math

# Step 1: Define the parameters based on the problem description.
# The problem uses a variable 'r' for the sidelength.
# The calculations are symbolic, so we can illustrate with r=1,
# as r will cancel out in the final ratio.
r_val = 1
num_arcs = 6
# The angle for each arc is pi/3 radians (60 degrees).
arc_angle_rad = math.pi / 3

# Step 2: Calculate the length of the curve C.
# The total length is the number of arcs times the length of one arc.
# Length of one arc = r * angle_in_radians
# Total length C = num_arcs * r * arc_angle_rad = 6 * r * (pi/3) = 2 * pi * r.
total_length_C_symbolic = "2 * pi * r"

# Step 3: Compute the required ratio.
# The problem asks to divide the length of C by 2 * pi * r.
# Ratio = (2 * pi * r) / (2 * pi * r)
ratio = (2 * math.pi * r_val) / (2 * math.pi * r_val)

# Step 4: Convert the ratio to a whole number percentage.
percentage = round(ratio * 100)

# Step 5: Print the results, showing the numbers in the final equation.
print(f"The curve C is composed of {num_arcs} arcs.")
print("Each arc has a length of (pi * r) / 3.")
print(f"The total length of C is given by the equation: Length(C) = {num_arcs} * (pi * r) / 3 = {total_length_C_symbolic}")
print("")
print("We need to compute the ratio: Length(C) / (2 * pi * r)")
print(f"Ratio = ({total_length_C_symbolic}) / (2 * pi * r) = {ratio:.0f}")
print("")
print(f"As a whole number percentage, the final answer is: {percentage}%")

# The final answer in the required format
# It is an integer, so we don't need a decimal point.
final_answer = int(percentage)
# print(f'<<<{final_answer}>>>')