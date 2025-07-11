import math

# Step 1: Define the given values based on the problem statement.
# S is the length of the string "SNYWR".
S_string = "SNYWR"
S = len(S_string)

# Side BC = a
a = S
# Angle B in degrees
angle_B_degrees = 30
# Convert angle B to radians for use in math functions
angle_B_radians = math.radians(angle_B_degrees)

# Step 2: Establish the condition for the ambiguous case (SSA) to have two solutions.
# For two possible triangles, the length of side AC (let's call it b) must be
# greater than the altitude from vertex C but less than the length of side BC.
# The condition is: a * sin(B) < b < a

# Step 3: Calculate the lower and upper bounds for the length of AC.
lower_bound = a * math.sin(angle_B_radians)
upper_bound = a

# Step 4: Display the logic and the resulting inequality with numerical values.
print(f"Given BC = a = {a} and angle B = {angle_B_degrees} degrees.")
print("For there to be two possible values for side AB, the length of AC must satisfy the inequality:")
print("BC * sin(B) < AC < BC")
# The following line shows the "final equation" with each number as requested.
print(f"{a} * sin({angle_B_degrees}) < AC < {a}")
print(f"{a * math.sin(angle_B_radians):.1f} < AC < {upper_bound}")

# Step 5: Find all integers within the calculated range for AC.
# The first integer greater than the lower bound
start_integer = math.floor(lower_bound) + 1
# The last integer less than the upper bound
end_integer = upper_bound - 1

integer_values_for_AC = []
if start_integer <= end_integer:
    for b in range(start_integer, end_integer + 1):
        integer_values_for_AC.append(b)

# Step 6: Print the list of possible integer values and the final count.
count = len(integer_values_for_AC)
print(f"\nThe possible integer values for AC are: {integer_values_for_AC}")
print(f"The number of integer values of AC for which there are two possible triangles is: {count}")