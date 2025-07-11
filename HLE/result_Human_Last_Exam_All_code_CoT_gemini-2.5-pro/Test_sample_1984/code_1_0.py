import math

# Step 1: Define the function f(x) from the differential equation x'(t) = f(x).
# f(x) = -x^3 + 2x^2 - x
# To find equilibrium points, we solve f(x) = 0.
# f(x) = -x(x^2 - 2x + 1) = -x(x - 1)^2
# The equilibrium points are x = 0 and x = 1.

# Step 2: Determine the stability of each point.
# A point x_e is stable if f(x) > 0 for x < x_e and f(x) < 0 for x > x_e.
#
# For x = 0:
# - Let's test x = -0.1 (a value less than 0): f(-0.1) = -(-0.1)(-0.1 - 1)^2 = 0.1 * (-1.1)^2 > 0.
# - Let's test x = 0.1 (a value greater than 0): f(0.1) = -(0.1)(0.1 - 1)^2 = -0.1 * (-0.9)^2 < 0.
# Since the sign changes from positive to negative, x = 0 is a stable equilibrium point.
#
# For x = 1:
# - Let's test x = 0.9 (a value less than 1): f(0.9) = -(0.9)(0.9 - 1)^2 = -0.9 * (-0.1)^2 < 0.
# - Let's test x = 1.1 (a value greater than 1): f(1.1) = -(1.1)(1.1 - 1)^2 = -1.1 * (0.1)^2 < 0.
# Since the function is negative on both sides of x = 1 (for values close to 1),
# the flow is away from x=1 on the left. Thus, x = 1 is not a stable equilibrium.

# Step 3: Count the number of stable equilibrium points.
# There is only one stable equilibrium point, so m = 1.
m = 1

# Step 4: Calculate the final expression m - 2^4048.
exponent = 4048
base = 2
second_term = base**exponent
result = m - second_term

# Step 5: Print the full equation with the calculated numbers.
print(f"{m} - {base}^{exponent} = {result}")