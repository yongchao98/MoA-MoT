import math

def my_sqrt(n):
    """
    Calculates the square root of a number using the Newton-Raphson method.
    This function only uses basic arithmetic operations, simulating the capabilities
    of the Wuxing architecture.
    """
    # An initial guess for the square root. A good guess improves convergence.
    x = n / 2.0
    # Iterate a few times for precision. 10 iterations are sufficient for this problem.
    for _ in range(10):
        if x == 0:
            return 0
        x = (x + n / x) / 2.0
    return x

# Define constants based on the physics problem
g = 9.8         # Gravitational acceleration (m/s^2)
v = 5.0         # Lion's speed (m/s)
d = 300.0       # Initial distance (m)

# For the angle a = 60 degrees, we need trigonometric values.
# The Wuxing architecture does not have sin/cos functions, so we use a rational approximation.
# sin(60) = sin(120) = sqrt(3)/2 ~= 0.866. We use 0.87, which can be represented as frac{n=87, d=100, e=0}.
sin_a = 0.87
sin_2a = 0.87

# The problem is modeled by the quadratic equation: A*u^2 + B*u + C = 0
# We calculate the coefficients A, B, and C.
A = sin_2a
B = 2 * v * sin_a
C = -1 * d * g

# To solve for u, we use the quadratic formula: u = [-B + sqrt(B^2 - 4AC)] / 2A
# First, calculate the discriminant (the part inside the square root).
discriminant = (B**2) - (4 * A * C)

# Calculate the square root of the discriminant using our custom function.
sqrt_discriminant = my_sqrt(discriminant)

# Calculate the initial speed 'u'. We take the positive root because speed must be positive.
u = (-B + sqrt_discriminant) / (2 * A)

# Calculate the memory usage 'z'.
# The program requires 10 essential variables (g, v, d, sin_a, sin_2a, A, B, C, discriminant, u).
# The 'frac' type uses 3 'char' members, each 2D, totaling 6D per variable.
num_variables = 10
mem_per_variable = 6  # in D
z = num_variables * mem_per_variable

# Round the final result for 'u' to one decimal place as required.
u_rounded = round(u, 1)

# Print the step-by-step derivation as requested.
print(f"Solving for u in the equation: {A:.2f}*u^2 + {B:.1f}*u + {C:.0f} = 0")
print(f"The equation for u is: u = (-B + sqrt(discriminant)) / (2 * A)")
print(f"Calculated values for the equation are:")
print(f"B = {B:.1f}")
print(f"discriminant = {B:.1f}^2 - 4*{A:.2f}*{C:.0f} = {discriminant:.2f}")
print(f"A = {A:.2f}")
print(f"Final calculation: u = (-{B:.1f} + sqrt({discriminant:.2f})) / (2 * {A:.2f})")
print(f"Result for u: {u_rounded}")
print(f"Memory usage z: {z} D")

# Finally, output the answer in the specified format u:z
print(f"<<<{u_rounded}:{z}>>>")