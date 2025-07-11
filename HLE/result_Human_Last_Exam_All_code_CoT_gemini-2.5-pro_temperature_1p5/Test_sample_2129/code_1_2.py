import math

# Parameters determined from the problem analysis
N = 17
lamb = 2
a = 5
pi = math.pi

# Derived value for y3(x0)
y3_x0 = (pi**4.5) / 320

# Final expression to calculate
final_value = (N + lamb) * (y3_x0)**(lamb / a)

# As requested, output each number in the final equation.
# The numbers are: (N + lamb), y3_x0, lamb, a
# (17 + 2) * ( (pi**4.5)/320 )**(2/5)
print(f"The equation to be solved is: ({N} + {lamb}) * ({y3_x0:.6f})**({lamb}/{a})")

# Print the result of the calculation
print(f"The result is: {final_value:.6f}")