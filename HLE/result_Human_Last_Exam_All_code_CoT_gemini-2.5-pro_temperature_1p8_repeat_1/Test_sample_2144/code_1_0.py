import math

# Step 1: The value of p at the target point y=-3.
# From solving the quartic equation p^4 - 18p - 27 = 0, we found an integer root.
p = 3

# We can verify that p=3 gives y=-3 with our parametric solution
# y(p) = -p^2 / sqrt(2p+3)
# y(3) = -(3^2) / sqrt(2*3+3) = -9 / sqrt(9) = -9/3 = -3. This confirms p=3 is the correct value.

# Step 2: Calculate x0 using the parametric equation for x.
# x(p) = -3 * (p + 1) / sqrt(2p + 3)
p_plus_1 = p + 1
two_p_plus_3 = 2 * p + 3
sqrt_two_p_plus_3 = math.sqrt(two_p_plus_3)
numerator = -3 * p_plus_1
x0 = numerator / sqrt_two_p_plus_3

# Step 3: Print the result, showing each number in the final equation.
print(f"The value of the derivative p=dy/dx at the target point is: {p}")
print(f"The corresponding position x0 is calculated using the formula: x0 = -3 * (p + 1) / sqrt(2*p + 3)")
print(f"Substituting p = {p}:")
print(f"x0 = -3 * ({p} + 1) / sqrt(2*{p} + 3)")
print(f"x0 = -3 * ({p_plus_1}) / sqrt({two_p_plus_3})")
print(f"x0 = {numerator} / {sqrt_two_p_plus_3}")
print(f"The position x0 is: {x0}")

# Final check of the point (x0, y=-3) with p=3 in the original ODE: (dy/dx)^3 + y^2 = xy(dy/dx)
lhs = p**3 + (-3)**2
rhs = x0 * (-3) * p
print(f"\nVerification:")
print(f"Checking if the point (x, y) = ({x0}, -3) with slope p = {p} satisfies the ODE:")
print(f"LHS = p^3 + y^2 = {p}^3 + (-3)^2 = {p**3} + {(-3)**2} = {lhs}")
print(f"RHS = x*y*p = {x0} * (-3) * {p} = {rhs}")
print(f"Since LHS = RHS ({lhs} = {rhs}), the solution is correct.")