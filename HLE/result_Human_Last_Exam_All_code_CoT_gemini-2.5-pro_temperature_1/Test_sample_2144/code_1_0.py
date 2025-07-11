import numpy as np
import math

# Step 1: Define the quartic equation for the parameter A.
# As derived from the problem's physics and substitutions, the parameter 'A'
# that defines the state of the particle must satisfy the equation:
# A^4 + 144*A - 432 = 0
# when the particle is at y = -3.
# We define the polynomial's coefficients to find the roots.
coeffs = [1, 0, 0, 144, -432]

# Step 2: Solve for the parameter A.
# We use numpy's root finding function to solve the quartic equation.
roots = np.roots(coeffs)

# The analysis shows that we should choose the integer root for the solution.
# We iterate through the roots to find the one that is a real integer.
A_val = None
for r in roots:
    # Check if the root is real and is close to an integer
    if np.isreal(r) and np.isclose(r.real, round(r.real)):
        A_val = round(r.real)
        break

# The integer root found is A = -6.

# Step 3: Calculate the position x0 from the parameter A.
# The x-coordinate is related to the parameter A by the formula:
# x = 3*(A-2) / (2 * sqrt(3-A))
# We substitute A = -6 into this formula.

numerator = 3 * (A_val - 2)
denominator = 2 * math.sqrt(3 - A_val)
x0 = numerator / denominator

# Step 4: Display the result and the final equation as requested.
# The final position x0 is where y(x0) = -3.
y_final = -3
x0_final = int(round(x0))

print(f"The position x0 where the particle reaches y = {y_final} has been calculated.")
print(f"The final state is described by the equation: y({x0_final}) = {y_final}")
print("Here are the numerical components of that final equation:")
print(f"The value of the final y-coordinate is: {y_final}")
print(f"The value of the final x-coordinate (x0) is: {x0_final}")
