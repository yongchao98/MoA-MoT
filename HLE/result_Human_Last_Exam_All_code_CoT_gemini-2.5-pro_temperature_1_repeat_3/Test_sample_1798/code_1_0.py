import numpy as np
from scipy.integrate import quad

# Step 1 & 2: Define the parameters and the function for the circle's upper arc.
# The circle passes through A(1,0), B(10,0), C(4,3).
# From the points, we derive the center (h, k) and radius squared (r2).
h = 5.5
k = -1.5
r2 = 22.5

def f(x):
    """
    Equation y=f(x) for the upper arc of the circle passing through A, B, and C.
    The equation is y = k + sqrt(r^2 - (x-h)^2).
    """
    # Ensure the value inside the square root is non-negative to avoid domain errors with quad
    radicand = r2 - (x - h)**2
    if isinstance(radicand, (np.ndarray, list, tuple)):
        radicand[radicand < 0] = 0
    elif radicand < 0:
        radicand = 0
    return k + np.sqrt(radicand)

# Step 3: Calculate alpha
# The total area under f(x) from x=1 to x=10 is the integral of f(x) over this domain.
# alpha = 1 / total_area
total_area, _ = quad(f, 1, 10)
alpha = 1 / total_area

# Step 4: Calculate P(X < 3)
# This is alpha * (integral of f(x) from 1 to 3).
area_lt_3, _ = quad(f, 1, 3)
prob_lt_3 = alpha * area_lt_3

# Final output
print("The equation of the circle is (x - 5.5)^2 + (y + 1.5)^2 = 22.5")
print(f"The value of alpha is: {alpha}")
print(f"The value of P(X < 3) is: {prob_lt_3}")