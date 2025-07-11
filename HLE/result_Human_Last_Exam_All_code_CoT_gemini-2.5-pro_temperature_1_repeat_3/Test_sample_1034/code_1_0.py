import math

# Step 1: Define the given values.
# S corresponds to SNYWR on a phone keypad.
S = 65294
angle_B_degrees = 30

# The side BC is 'a' in the standard a, b, c notation.
a = S
# The side AC is 'b'.

# Step 2: Calculate the bounds for AC based on the ambiguous case condition.
# The condition for two possible triangles is a * sin(B) < b < a.
# sin(30 degrees) is exactly 0.5.
lower_bound = a * 0.5
upper_bound = a

print(f"The given side BC = {a}.")
print(f"The given angle B = {angle_B_degrees} degrees.")
print("For two possible triangles to exist, the length of side AC must be in the range (BC * sin(B), BC).")
print(f"This gives the inequality: {lower_bound} < AC < {upper_bound}")
print("-" * 30)

# Step 3: Count the number of integers in this range.
# The first integer value for AC is the integer part of the lower bound + 1.
first_integer_AC = int(lower_bound) + 1

# The last integer value for AC is the upper bound - 1.
last_integer_AC = upper_bound - 1

print(f"The smallest possible integer value for AC is {first_integer_AC}.")
print(f"The largest possible integer value for AC is {last_integer_AC}.")

# Step 4: Calculate the total count.
# The number of integers is (last - first + 1).
# The final equation to calculate the count is:
print(f"Final count equation: {last_integer_AC} - {first_integer_AC} + 1")
number_of_values = last_integer_AC - first_integer_AC + 1

print(f"The number of integer values of AC is: {number_of_values}")