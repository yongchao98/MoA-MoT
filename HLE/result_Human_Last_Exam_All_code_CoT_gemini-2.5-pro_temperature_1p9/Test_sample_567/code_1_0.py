import math

# The critical value of 'a' is determined by the fourth power of the golden ratio, phi.
# The formula for the golden ratio is: phi = (1 + sqrt(5)) / 2
# Therefore, the equation for 'a' is: a = ((1 + sqrt(5)) / 2)^4

# Define the numbers that constitute the final equation.
num_one = 1
num_five = 5
num_two = 2
power = 4

# Calculate the golden ratio, phi.
phi = (num_one + math.sqrt(num_five)) / num_two

# Calculate 'a' by raising phi to the specified power.
a = phi**power

# As requested, we will now output each number in the final equation
# by printing the full equation with its components and the final result.
print(f"The equation for the critical value of 'a' is: a = (({num_one} + sqrt({num_five})) / {num_two}) ^ {power}")
print(f"The calculated value is: a = {a}")
