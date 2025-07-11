import math

# Define the constants for the problem
n = 2
u1 = 0.495  # Lower bound of the interval (x1/a)
u2 = 0.505  # Upper bound of the interval (x2/a)

# The analytical solution for the probability P is:
# P = [u - (1/(2*n*pi)) * sin(2*n*pi*u)] evaluated from u1 to u2
# P = (u2 - u1) - (1/(2*n*pi)) * (sin(2*n*pi*u2) - sin(2*n*pi*u1))

# Calculate the arguments for the sine functions
arg1 = 2 * n * math.pi * u1
arg2 = 2 * n * math.pi * u2

# Calculate the sine terms
sin_term1 = math.sin(arg1)
sin_term2 = math.sin(arg2)

# Calculate the different parts of the equation
delta_u = u2 - u1
prefactor = 1 / (2 * n * math.pi)

# Calculate the final probability
probability = delta_u - prefactor * (sin_term2 - sin_term1)

# Print the equation with the numbers plugged in
print("The probability P is calculated from the equation:")
print(f"P = (u2 - u1) - (1/(2*n*pi)) * [sin(2*n*pi*u2) - sin(2*n*pi*u1)]")
print("\nPlugging in the values:")
print(f"P = ({u2} - {u1}) - (1/(2*{n}*pi)) * [sin(2*{n}*pi*{u2}) - sin(2*{n}*pi*{u1})]")
print(f"P = {delta_u} - ({prefactor:.6f}) * [{sin_term2:.6f} - ({sin_term1:.6f})]")
print("\nFinal Result:")
print(f"The probability is: {probability}")
