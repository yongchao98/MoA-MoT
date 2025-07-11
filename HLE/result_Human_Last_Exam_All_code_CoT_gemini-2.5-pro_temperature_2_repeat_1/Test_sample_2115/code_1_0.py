import math

# This script calculates the value of the definite integral.

# Value of the denominator at x=1: e^2 + e + 1
val_at_1 = math.e**2 + math.e + 1

# Value of the denominator at x=0: e^0 + e^0 + 1 = 1 + 1 + 1 = 3
val_at_0 = 3.0

# The analytical solution to the integral is 3 * ln(3 / (e^2 + e + 1))
result = 3 * math.log(val_at_0 / val_at_1)

# Print out the equation with the calculated values
print("The final equation is: 3 * (ln(3) - ln(e^2 + e + 1))")
print(f"Which evaluates to: 3 * (ln({val_at_0}) - ln({math.e**2} + {math.e} + 1))")
print(f"                   = 3 * ({math.log(val_at_0)} - {math.log(val_at_1)})")
print(f"                   = 3 * ({math.log(val_at_0) - math.log(val_at_1)})")
print(f"                   = {3 * (math.log(val_at_0) - math.log(val_at_1))}")

# Print the final numerical answer
print("\nThe numerical result of the integral is:")
print(result)