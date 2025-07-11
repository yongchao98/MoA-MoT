import math

# Step 1: Define the value of S based on the SNYWR placeholder.
# Using a phone keypad mapping: S=7, N=6, Y=9, W=9, R=7
S = 76997

# Step 2: Define the given side and angle.
# Side BC is given as S.
# Angle B is 30 degrees, so sin(B) = 0.5.
BC = S
angle_B_degrees = 30

# Step 3: The condition for two possible triangles in the SSA case is
# BC * sin(B) < AC < BC.
# This translates to S * 0.5 < AC < S.
lower_bound = S / 2
upper_bound = S

# Step 4: Find the number of integers for AC in the range (lower_bound, upper_bound).
# The smallest integer AC can be is the first integer greater than the lower bound.
first_integer_AC = math.floor(lower_bound) + 1

# The largest integer AC can be is the last integer less than the upper bound.
last_integer_AC = math.ceil(upper_bound) - 1

# Step 5: Calculate the total count of possible integer values.
# The number of integers in a range [first, last] is (last - first + 1).
count = last_integer_AC - first_integer_AC + 1

# Print the explanation and the final equation as requested.
print(f"The given side length is S = {S}.")
print(f"For two triangles to exist, the length of side AC must be in the range ({lower_bound}, {upper_bound}).")
print("We need to count the number of integers in this range.")
print(f"The smallest integer value for AC is {first_integer_AC}.")
print(f"The largest integer value for AC is {last_integer_AC}.")
print("\nThe final equation to find the number of integer values is (Last - First + 1):")
print(f"{last_integer_AC} - {first_integer_AC} + 1 = {count}")