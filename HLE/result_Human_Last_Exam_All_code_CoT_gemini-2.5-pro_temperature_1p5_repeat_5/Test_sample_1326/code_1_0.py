import numpy as np
from scipy.integrate import quad
import warnings

def p1(x):
    """
    The first piece of the function, defined for 0 <= x <= 3.
    """
    return (2 * x**3) / 8

def p2(x):
    """
    The second piece of the function, defined for 3 <= x <= 5.
    """
    return np.exp(x) * (1 + np.sin(x)) / (1 + np.cos(x))

# The integration range is [0, 4]. We split it at x=3.

# 1. Integrate the first part from 0 to 3
integral1, error1 = quad(p1, 0, 3)

# 2. Integrate the second part from 3 to 4.
# This part has a singularity at x=pi, which is in the interval [3, 4].
# We inform quad about this singularity using the 'points' argument.
# This will calculate the Cauchy Principal Value.
# We also catch and print any warnings from the integration routine.
print("Calculating the integral for the second part of the function from x=3 to x=4.")
print("Note: The integrand has a singularity at x=pi (approx 3.14159) within this interval.")

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    
    integral2, error2 = quad(p2, 3, 4, points=[np.pi])
    
    if w:
        print("\n--- Warnings from SciPy ---")
        for warning_message in w:
            print(f"Warning: {warning_message.message}")
        print("---------------------------\n")

# 3. Sum the results
total_integral = integral1 + integral2

# Print the final result in the format of an equation
print("The value of the integral from 0 to 3 is:", integral1)
print("The principal value of the integral from 3 to 4 is:", integral2)
print("The total integral is the sum of these two values:")
print(f"{integral1} + ({integral2}) = {total_integral}")

# Final numeric answer
print("\nFinal Answer:")
print(total_integral)