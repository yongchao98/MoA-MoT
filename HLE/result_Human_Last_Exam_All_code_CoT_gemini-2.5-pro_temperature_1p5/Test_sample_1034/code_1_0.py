import math

# Step 1: Assign the given values.
# S = SNYWR maps to 26734.
BC = 26734
angle_B_degrees = 30

# Step 2: Define the condition for two possible triangles.
# In triangle ABC, let a = BC, b = AC. The condition for two solutions is h < b < a,
# where h is the altitude from C, given by h = a * sin(B).

# Step 3: Calculate the bounds for the length of AC.
# Convert angle B to radians for the sin function.
angle_B_radians = math.radians(angle_B_degrees)
# a = BC
a = BC
# Calculate the altitude h, which is the lower bound for AC.
h = a * math.sin(angle_B_radians)

# The bounds for AC are (h, a).
lower_bound = h
upper_bound = a

print(f"For two possible triangles, the length of AC must be between the altitude h and the length of side BC.")
print(f"Side BC = {upper_bound}")
print(f"Altitude h = BC * sin(B) = {a} * sin({angle_B_degrees}°) = {lower_bound}")
print(f"So, the condition is {lower_bound} < AC < {upper_bound}.")
print("")

# Step 4: Count the number of integer values for AC.
# The smallest integer value for AC is floor(h) + 1.
# The largest integer value for AC is ceil(a) - 1.
first_integer_AC = int(lower_bound) + 1
last_integer_AC = int(upper_bound) - 1

# Calculate the total number of integers.
count = last_integer_AC - first_integer_AC + 1

print(f"The number of possible integer values for AC is the count of integers in the range ({first_integer_AC}, {last_integer_AC}), inclusive.")
print("The calculation is:")
# Output the final equation as requested
print(f"Number of values = {last_integer_AC} - {first_integer_AC} + 1 = {count}")
