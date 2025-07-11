import math

# We are given the system of differential equations:
# x'(t) = -3*x(t)*y(t)
# y'(t) = -y(t)^2 - x(t) + 1
# with the initial condition x(0) > 1.

# Through phase plane analysis, we find that the behavior of the solution
# is separated by a curve (separatrix) that passes through the saddle
# point (1, 0). The equation of this separatrix is:
# y^2 = 2*x + 1 - 3*x^(2/3)

# For a solution to blow up (i.e., y(t) -> -infinity in finite time),
# the initial point (x(0), y(0)) must lie below the stable manifold of
# this separatrix. This gives the condition:
# y(0) < sqrt(2*x(0) + 1 - 3*x(0)^(2/3))

# Let's demonstrate this with a specific value for x(0). We'll choose
# x(0) = 8, which satisfies the condition x(0) > 1 and has a nice integer cube root.
x0 = 8

# The parameters of the critical value equation are:
c1 = 2
c2 = 1
c3 = 3
power_numerator = 2
power_denominator = 3
power = power_numerator / power_denominator

# Calculate the term inside the square root
val_inside_sqrt = c1 * x0 + c2 - c3 * (x0**power)

# Calculate the critical value for y(0)
critical_y0 = math.sqrt(val_inside_sqrt)

# As requested, we will output each number in the final equation.
print(f"To find the range of y(0) for blow-up, we use the inequality:")
print(f"y(0) < sqrt(c1 * x(0) + c2 - c3 * x(0)^(num/den))")
print(f"The numbers in the equation are:")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print(f"c3 = {c3}")
print(f"num = {power_numerator}")
print(f"den = {power_denominator}")

print(f"\nFor the chosen initial condition x(0) = {x0}, we substitute the values:")
print(f"y(0) < sqrt({c1} * {x0} + {c2} - {c3} * {x0}^({power_numerator}/{power_denominator}))")
print(f"y(0) < sqrt({c1 * x0} + {c2} - {c3} * {x0**power})")
print(f"y(0) < sqrt({val_inside_sqrt})")

print("\nTherefore, for x(0) = 8, the solution blows up if y(0) satisfies:")
print(f"y(0) < {critical_y0}")
