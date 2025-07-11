import math

# --- Step 1: Define constants of the problem ---
d = 300.0  # Initial distance in meters
v = 5.0    # Lion's speed in m/s
a_deg = 60.0 # Angle of throw in degrees
g = 9.8    # Acceleration due to gravity in m/s^2

# Convert angle to radians for use in math functions
a_rad = math.radians(a_deg)

# --- Step 2: Formulate the quadratic equation Au^2 + Bu + C = 0 ---
# The equation for the initial speed 'u' is:
# u^2 * sin(2a) + u * (2*v*sin(a)) - d*g = 0

# Calculate the coefficients A, B, and C
A = math.sin(2 * a_rad)
B = 2 * v * math.sin(a_rad)
C = -d * g

# --- Step 3: Print the equation with its numerical coefficients ---
print("The problem is solved using the quadratic equation: A*u^2 + B*u + C = 0")
print("Where 'u' is the initial speed of the rock.")
print("\nCalculating the coefficients based on the given values:")
print(f"A = sin(2 * {a_deg}°) = {A:.4f}")
print(f"B = 2 * {v} * sin({a_deg}°) = {B:.4f}")
print(f"C = -{d} * {g} = {C:.4f}")
print("\nSo, the final equation is:")
# The formatting ensures that the sign of B and C is displayed correctly.
print(f"{A:.4f}u^2 + {B:.4f}u {C:.4f} = 0")


# --- Step 4: Solve the quadratic equation for u ---
# Using the quadratic formula: u = (-B ± sqrt(B^2 - 4AC)) / 2A
# We take the positive root because speed must be positive.
discriminant = B**2 - 4 * A * C
u = (-B + math.sqrt(discriminant)) / (2 * A)

# --- Step 5: Calculate memory usage 'z' ---
# Variables needed: d, v, g, sin(a), cos(a), u -> 6 variables
# Data type: frac (size = 2D + 2D + 2D = 6D)
# Total memory z = 6 variables * 6D/variable
z = 6 * 6

# --- Step 6: Display the final answer ---
print("\nSolving for u gives the initial speed:")
print(f"u = {u:.1f} m/s")

print("\nCalculating memory usage 'z':")
print(f"z = {z} D")

print("\nFinal Answer (u:z):")
print(f"{u:.1f}:{z}")

# <<<53.5:36>>>