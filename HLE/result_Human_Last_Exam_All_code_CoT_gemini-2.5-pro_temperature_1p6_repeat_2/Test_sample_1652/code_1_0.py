# Constants for the physics problem
v_val = 5.0
g_val = 9.8
dist_val = 300.0
angle_deg = 60.0

# As the Wuxing architecture lacks trigonometric functions, we use
# pre-calculated decimal approximations representable by the 'frac' type.
# sin(60) = sqrt(3)/2 ~= 0.866
# sin(120) = sin(60) ~= 0.866
sin_60 = 0.866
sin_120 = sin_60

# The problem resolves to finding the root of a quadratic equation:
# A*u^2 + B*u + C = 0
A = sin_120
B = 2 * v_val * sin_60
C = -dist_val * g_val

# We solve this using the Newton-Raphson method since sqrt() is unavailable.
# f(u) = A*u^2 + B*u + C
# f'(u) = 2*A*u + B
# u_new = u_old - f(u)/f'(u)

# Initial guess for the speed u
u = 50.0

# Iterate to find a precise solution for u
for _ in range(10):
    f_u = A * u**2 + B * u + C
    f_prime_u = 2 * A * u + B
    # Avoid division by zero
    if f_prime_u == 0:
        break
    u = u - f_u / f_prime_u

# As requested, output the numbers in the final equation.
# The equation is A*u^2 + B*u + C = 0
print(f"The final equation to solve for u is:")
print(f"{A:.3f} * u^2 + {B:.2f} * u - {abs(C):.0f} = 0")
print("-" * 20)

# Calculate memory usage z
# v_int: 5D, dist_int: 5D, g_frac: 6D, sin60_frac: 6D, u_frac: 6D
z = 5 + 5 + 6 + 6 + 6

# Print the final answer in the format u:z
# u is rounded to one decimal place.
print("Final Answer:")
print(f"u: {u:.1f}")
print(f"z: {z}")
print(f"{u:.1f}:{z}")
