import math

# Step 1 & 2: Define problem parameters based on the Wuxing architecture constraints.
# All values that would be `frac` types in the Wuxing C compiler are represented here.

# v = 5 m/s
v = 5.0
# Initial distance = 300 m
dist = 300.0
# g = 9.8 m/s^2. We use floats here for the Python simulation.
g = 9.8
# Angle a = 60 degrees.
# Wuxing C does not have a sin function. We must use a fractional approximation.
# sin(60) = sin(120) ≈ 0.866. A good fractional approximation is 7/8 = 0.875.
sin_a = 7.0 / 8.0
sin_2a = 7.0 / 8.0 # sin(120) = sin(180-60) = sin(60)

# Step 3: Use an iterative approach to find u.
# We are solving for u in the equation: u^2*sin(2a) + (2*v*u*sin(a)) - dist*g = 0
# Or, checking when: u^2*sin(2a) + 2*v*u*sin_a >= dist*g

rhs = dist * g
u = 0.0
step = 0.1 # Iterate with a precision of 0.1

while True:
    u += step
    lhs = (u * u * sin_2a) + (2 * v * u * sin_a)

    if lhs >= rhs:
        u_final = round(u, 1)
        break

# Step 4: Calculate the memory usage 'z'.
# The program requires storage for the main variables. In the Wuxing architecture,
# these would be of type 'frac'.
# Each 'frac' variable (n, d, e) costs 2D + 2D + 2D = 6D.
# The variables are: g, v, dist, sin_a, u, step.
# Total variables = 6
num_variables = 6
cost_per_variable_D = 6
z = num_variables * cost_per_variable_D

# Step 5: Output the result in the format u:z.
# The problem asks to output the numbers in the final equation.
# Our "final equation" is the required answer format u:z.
print(f"{u_final}:{z}")
<<<53.2:36>>>