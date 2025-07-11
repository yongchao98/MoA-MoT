import numpy as np
from scipy.integrate import quad

# Step 1: Define the parameters of the circle equation
# (x - h)^2 + (y - k)^2 = r2
h = 5.5
k = -1.5
r2 = 22.5

# Step 2: Define the function f(x)
def f(x):
    """
    Equation of the upper arc of the circle passing through A, B, and C.
    """
    # Ensure the term inside the square root is non-negative
    # This is guaranteed for x in [1, 10]
    discriminant = r2 - (x - h)**2
    return k * -1 + np.sqrt(discriminant)

# Step 3: Calculate the normalization constant alpha
# alpha = 1 / integral of f(x) from 1 to 10
total_integral, _ = quad(f, 1, 10)
alpha = 1 / total_integral

# Step 4: Calculate the probability P(X < 3)
# P(X < 3) = alpha * integral of f(x) from 1 to 3
partial_integral, _ = quad(f, 1, 3)
probability_x_lt_3 = alpha * partial_integral

# Print the final results
print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2")
print(f"where h = {h}, k = {k}, and r^2 = {r2}")
print("This can be written as:")
print(f"(x - {h})^2 + (y - ({k}))^2 = {r2}")
print("\n---")
print(f"The value of alpha is: {alpha}")
print(f"The value of P(X < 3) is: {probability_x_lt_3}")
