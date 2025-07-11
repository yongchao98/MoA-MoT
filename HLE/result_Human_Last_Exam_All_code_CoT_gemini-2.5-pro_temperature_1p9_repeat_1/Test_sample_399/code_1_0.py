# The explicit solution to the initial value problem is x(t) = -8 / (4*t^2 - 8*t + 1).
# We need to find the value of x(1).

# Set the value of t
t = 1

# Define the terms in the formula
numerator = -8
# We will calculate the denominator: 4*t^2 - 8*t + 1
term1 = 4 * (t**2)
term2 = -8 * t
term3 = 1
denominator = term1 + term2 + term3

# Calculate the final result
result = numerator / denominator

# Print the final equation with the numbers plugged in
print("The solution formula is x(t) = -8 / (4*t^2 - 8*t + 1)")
print("Substituting t=1:")
print(f"x(1) = {numerator} / (4 * ({t})^2 - 8 * ({t}) + 1)")
print(f"x(1) = {numerator} / ({term1} {term2} + {term3})")
print(f"x(1) = {numerator} / {denominator}")
print(f"x(1) = {result}")
