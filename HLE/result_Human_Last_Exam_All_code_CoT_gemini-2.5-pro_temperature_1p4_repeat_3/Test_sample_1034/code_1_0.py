import math

# Step 1 & 4: Define the given values.
# S is derived from SNYWR using a phone keypad mapping (S=7, N=6, Y=9, W=9, R=7).
S = 76997
# Given angle B in degrees.
angle_B_degrees = 30

# Step 2 & 3: Formulate the inequality for the ambiguous case.
# The condition for two possible triangles is S * sin(B) < AC < S.
# We need to find the number of integers for AC in this range.

# Step 5: Solve the inequality.
# sin(30 degrees) = 0.5
sin_B = 0.5
lower_bound = S * sin_B
upper_bound = S

# Step 6: Count the integers in the range (lower_bound, upper_bound).
# We need to find the number of integers AC such that lower_bound < AC < upper_bound.
# The smallest integer value for AC is the first integer greater than lower_bound.
first_integer_ac = math.floor(lower_bound) + 1

# The largest integer value for AC is the first integer less than upper_bound.
last_integer_ac = math.ceil(upper_bound) - 1

# The total number of integer values is (last - first + 1).
count = last_integer_ac - first_integer_ac + 1

# Output the results as requested.
print(f"The condition for two possible triangles is: {lower_bound} < AC < {upper_bound}")
print(f"The smallest integer value for AC is: {first_integer_ac}")
print(f"The largest integer value for AC is: {last_integer_ac}")
print("The number of integer values for AC is calculated as:")
print(f"{last_integer_ac} - {first_integer_ac} + 1 = {count}")