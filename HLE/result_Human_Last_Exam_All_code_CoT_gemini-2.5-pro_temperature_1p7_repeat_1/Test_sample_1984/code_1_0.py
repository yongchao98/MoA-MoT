import math

# Step 1: Find the number of stable equilibrium points, m.
# The differential equation is x'(t) = -x^3 + 2x^2 - x.
# To find equilibrium points, we set x'(t) = 0:
# -x^3 + 2x^2 - x = 0
# -x(x^2 - 2x + 1) = 0
# -x(x - 1)^2 = 0
# The equilibrium points are x = 0 and x = 1.

# Step 2: Determine the stability of these points.
# Let f(x) = -x^3 + 2x^2 - x. The stability is determined by the sign of f'(x).
# f'(x) = -3x^2 + 4x - 1.
# For x = 0: f'(0) = -1 < 0, so x=0 is a stable equilibrium point.
# For x = 1: f'(1) = -3(1)^2 + 4(1) - 1 = 0. The test is inconclusive.
# We must analyze the sign of f(x) near x=1.
# f(x) = -x(x-1)^2. For x near 1, x>0 and (x-1)^2 is always non-negative.
# Thus, f(x) is less than or equal to 0 near x=1.
# If x > 1, x'(t) < 0, so x(t) moves toward 1.
# If x < 1, x'(t) < 0, so x(t) moves away from 1.
# Since trajectories starting on one side of x=1 move away from it, x=1 is not a stable equilibrium point (it is semi-stable).
# Therefore, there is only one stable equilibrium point.
m = 1

# Step 3: Calculate the final expression m - 2^4048.
base = 2
exponent = 4048

# Python's integers have arbitrary precision, so it can handle this large number.
result = m - (base**exponent)

# The prompt asks to output each number in the final equation.
# So we print the equation with the values substituted.
print(f"{m} - {base}**{exponent} = {result}")