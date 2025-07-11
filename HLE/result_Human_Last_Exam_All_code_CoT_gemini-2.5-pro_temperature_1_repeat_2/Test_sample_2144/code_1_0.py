import math

# Based on the derivation, the specific solution to the differential equation
# satisfying the initial condition y(0) = -1 is given by the parametric equations:
# x(p) = -(3*p + 3) / sqrt(2*p + 3)
# y(p) = -p^2 / sqrt(2*p + 3)
# where p = dy/dx and p > -1.5.

# We need to find the position x0 where y(x0) = -3.
# We set y(p) = -3:
# -3 = -p^2 / sqrt(2*p + 3)
# This simplifies to the quartic equation: p^4 - 18*p - 27 = 0.

# By inspection or using a numerical solver, we can find the integer root of this equation.
# Let's test integer factors of 27. Testing p=3:
# 3^4 - 18*3 - 27 = 81 - 54 - 27 = 0.
# So, p = 3 is a valid solution for the parameter p.

p_final = 3

# Now, we calculate the corresponding x0 using the equation for x(p).
x0 = -(3 * p_final + 3) / math.sqrt(2 * p_final + 3)

print(f"The value of the parameter p where y = -3 is: {p_final}")
print(f"The corresponding position x0 is: {x0}")

# Finally, let's verify the solution by plugging the values into the original differential equation:
# (dy/dx)^3 + y^2 = x*y*(dy/dx)

y0 = -3
dy_dx = p_final

# Calculate the left-hand side (LHS) of the equation
lhs = dy_dx**3 + y0**2

# Calculate the right-hand side (RHS) of the equation
rhs = x0 * y0 * dy_dx

print("\nVerification of the solution:")
print(f"Plugging x = {x0}, y = {y0}, and dy/dx = {dy_dx} into the original equation:")
print(f"({dy_dx})^3 + ({y0})^2 = ({x0}) * ({y0}) * ({dy_dx})")
print(f"{dy_dx**3} + {y0**2} = {x0 * y0} * {dy_dx}")
print(f"{lhs} = {rhs}")

if math.isclose(lhs, rhs):
    print("The equation holds true. The solution is correct.")
else:
    print("The equation does not hold. There is an error in the solution.")
