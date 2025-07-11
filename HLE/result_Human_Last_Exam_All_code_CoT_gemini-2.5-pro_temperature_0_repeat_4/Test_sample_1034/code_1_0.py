import math

# Step 1: Define the problem variables
# Let S = SNYWR. We interpret this as the sum of the alphabetical positions.
s_str = "SNYWR"
# ord(c) - ord('A') + 1 gives the position in the alphabet (A=1, B=2, ...)
side_bc = sum(ord(c) - ord('A') + 1 for c in s_str)

angle_b_deg = 30
angle_b_rad = math.radians(angle_b_deg)

# Step 2: Establish the condition for two possible triangles (ambiguous SSA case)
# The condition is: side_bc * sin(B) < side_ac < side_bc
# This ensures that side AC is long enough to reach the base, but short enough to pivot to two positions.

# Step 3: Calculate the bounds for the length of side AC
lower_bound = side_bc * math.sin(angle_b_rad)
upper_bound = side_bc

# Step 4: Find the number of integers within this range
# The first possible integer for AC is the integer just above the lower bound.
first_integer_ac = math.floor(lower_bound) + 1

# The last possible integer for AC is the integer just below the upper bound.
last_integer_ac = math.ceil(upper_bound) - 1

# The total number of integers is the difference, plus one.
count = last_integer_ac - first_integer_ac + 1

# Step 5: Print the results as requested
print(f"Given BC = S = {side_bc} and angle B = {angle_b_deg} degrees.")
print("For two possible values of side AB, the length of side AC must satisfy the inequality:")
print(f"BC * sin(B) < AC < BC")
print("\nThe final equation with the calculated values is:")
# The prompt requires printing each number in the final equation.
print(f"{lower_bound} < AC < {upper_bound}")
print(f"\nThe range of integer values for AC is from {first_integer_ac} to {last_integer_ac}.")
print(f"\nThe number of integer values for AC is: {count}")
