import math

# We are given x(0) > 1. For this script, we'll use a specific value for x(0).
x0 = 2.0

# The critical value for y(0) is determined by the separatrix equation, which separates
# trajectories that blow up from those that do not.
# The blow-up condition is y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).

# We calculate the components of the expression under the square root.
term1 = 2 * x0
term2 = 1.0
term3 = 3 * (x0**(2.0/3.0))

# Calculate the value inside the square root.
val_inside_sqrt = term1 + term2 - term3

# For x0 > 1, the value inside the square root is positive, so a real solution exists.
# We can now calculate the critical threshold for y(0).
y_crit = math.sqrt(val_inside_sqrt)

# Display the reasoning and the result.
print(f"For the given system with the initial condition x(0) = {x0}, the solution blows up if y(0) is below a critical value.")
print("This critical value is derived from the separatrix equation y^2 = 2*x + 1 - 3*x^(2/3).")
print("\nWe evaluate the right-hand side at x = x(0):")
print(f"y_crit^2 = 2*({x0}) + {term2} - 3*({x0})^(2/3)")
print(f"y_crit^2 = {term1} + {term2} - {term3:.6f}")
print(f"y_crit^2 = {val_inside_sqrt:.6f}")
print(f"y_crit = {y_crit:.6f}")

print("\nThus, for the solution to blow up, the initial condition y(0) must satisfy:")
print(f"y(0) < {y_crit:.6f}")

print(f'<<<y(0) < {y_crit:.6f}>>>')