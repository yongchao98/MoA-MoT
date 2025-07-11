import math

# Step 1: Define the value for S from the problem statement.
# S = SNYWR is a placeholder. We will use a representative 5-digit number.
S = 54321

# Step 2: Establish the known values from the problem.
BC = S
angle_B_degrees = 30

print(f"Given S = {S}, side BC = {S}, and angle B = 30 degrees.")
print("For there to be two possible triangles, the length of side AC must be greater than the altitude from vertex C and less than the length of side BC.")
print("-" * 20)

# Step 3: Calculate the bounds for the length of side AC.
# The altitude h = BC * sin(B) = S * sin(30) = S / 2.
lower_bound = S / 2.0
upper_bound = S

print("The condition for AC is: altitude < AC < BC")
print(f"Calculating the altitude: h = {S} * sin(30) = {S} * 0.5 = {lower_bound}")
print(f"So, the inequality for AC is: {lower_bound} < AC < {upper_bound}")
print("-" * 20)

# Step 4: Find the number of integers for AC in this range.
# The smallest integer greater than the lower bound.
first_integer = math.floor(lower_bound) + 1

# The largest integer less than the upper bound.
last_integer = upper_bound - 1

# Calculate the total count of integers.
count = last_integer - first_integer + 1

print("We need to find the number of integer values for AC in this range.")
print(f"The smallest possible integer for AC is {first_integer}.")
print(f"The largest possible integer for AC is {last_integer}.")
print(f"The total number of integer values for AC is the count of integers from {first_integer} to {last_integer}.")
print(f"Calculation: {last_integer} - {first_integer} + 1 = {count}")
print("-" * 20)
print(f"Final Answer: The number of integer values for AC is {count}.")
