import math

# The value at which to evaluate the PDF
z = 0.2

# The PDF of Z is given by the formula f(z) = -2 * (z*ln(z) + (1-z)*ln(1-z)).
# We need to calculate f(0.2).

# Let's break down the calculation for the printout.
val_z = z
val_1_minus_z = 1 - z
log_z = math.log(val_z)
log_1_minus_z = math.log(val_1_minus_z)

# Calculate the final result
result = -2 * (val_z * log_z + val_1_minus_z * log_1_minus_z)

# Output the explanation and the step-by-step calculation
print("The problem simplifies to finding the PDF of a variable Z chosen uniformly from an interval [X, Y],")
print("where X and Y are independent random variables uniformly distributed on [0, 1].")
print("The resulting PDF is f(z) = -2 * (z*ln(z) + (1-z)*ln(1-z)).\n")
print(f"We want to calculate the value of this function at z = {z}.")
print("\n--- Calculation Steps ---")
print(f"f({z}) = -2 * ({z} * ln({z}) + {1-z} * ln({1-z}))")
print(f"f({z}) = -2 * ({z} * ({log_z:.5f}) + {1-z} * ({log_1_minus_z:.5f}))")
print(f"f({z}) = -2 * ({(val_z * log_z):.5f} + {(val_1_minus_z * log_1_minus_z):.5f})")
print(f"f({z}) = -2 * ({(val_z * log_z + val_1_minus_z * log_1_minus_z):.5f})")
print(f"f({z}) = {result:.5f}\n")
print("The exact analytical value is 2*ln(5) - 3.2*ln(2).")
print(f"The final numerical value of f(0.2) is: {result}")