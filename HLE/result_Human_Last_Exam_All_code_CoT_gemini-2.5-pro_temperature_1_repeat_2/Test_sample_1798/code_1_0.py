import numpy as np
from scipy.integrate import quad

# Step 1 & 2: Define the function f(x) based on the circle's equation.
# Circle center (h, k) and radius squared r_sq
h = 5.5
k = -1.5
r_sq = 22.5

def f(x):
    """Equation of the upper arc of the circle."""
    return k + np.sqrt(r_sq - (x - h)**2)

# Step 3: Find alpha
# The total integral is the normalization constant for the PDF.
total_integral, _ = quad(f, 1, 10)

# alpha is the reciprocal of the total integral.
alpha = 1 / total_integral

# Step 4: Calculate P(X < 3)
# The probability is the ratio of the partial integral to the total integral.
# The "final equation" is P(X<3) = partial_integral / total_integral
partial_integral, _ = quad(f, 1, 3)

probability = partial_integral / total_integral

# Output the results
print("Equation of the circle: (x - 5.5)^2 + (y + 1.5)^2 = 22.5")
print(f"The function is f(x) = -1.5 + sqrt(22.5 - (x - 5.5)^2)")
print("\n--- Calculating alpha ---")
print(f"The integral of f(x) from 1 to 10 is: {total_integral}")
print(f"alpha = 1 / (integral of f(x)) = {alpha}")

print("\n--- Calculating P(X < 3) ---")
print(f"The integral of f(x) from 1 to 3 is: {partial_integral}")
print(f"P(X < 3) = (integral from 1 to 3) / (integral from 1 to 10) = {probability}")

print("\nFinal values:")
print(f"alpha = {alpha}")
print(f"P(X < 3) = {probability}")