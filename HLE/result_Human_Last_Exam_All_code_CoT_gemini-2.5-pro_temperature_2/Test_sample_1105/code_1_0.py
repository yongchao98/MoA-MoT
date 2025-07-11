import math

# --- Step 1: Define physical constants and initial conditions ---
# Note: The problem requires a finite l_max, which means the spaceship must eventually stop.
# The parameters below, especially the low initial velocity v_0, are chosen to ensure this condition is met.

# Gravitational constant in N*(m/kg)^2
G = 6.674e-11
# Mass of the asteroid A (m) in kg (similar to a small moon)
m = 5.972e18
# Mass of the spaceship B (M) in kg
M = 1.0e5
# Initial distance between A and B (l_0) in meters
l_0 = 1.0e6
# Initial speed of the spaceship (v_0) in m/s
v_0 = 12.0
# Additional constant force on the spaceship (F) in Newtons
F = 10.0

# --- Step 2: Formulate and calculate the coefficients of the quadratic equation ---
# From the work-energy theorem, we derive an equation of the form:
# a*l_max^2 + b*l_max + c = 0

# Coefficient 'a'
a = F
# Coefficient 'b'
b = (0.5 * M * v_0**2) - (F * l_0) - (G * m * M / l_0)
# Coefficient 'c'
c = G * m * M

# --- Step 3: Print the breakdown of the equation and solve it ---
print("The maximum distance l_max can be found by solving a quadratic equation derived from the work-energy theorem.")
print(f"\nThe equation is of the form: a*(l_max)^2 + b*(l_max) + c = 0")
print("\nBased on the initial parameters:")
print(f"  G = {G:.3e} N(m/kg)^2, m = {m:.3e} kg, M = {M:.3e} kg")
print(f"  l_0 = {l_0:.3e} m, v_0 = {v_0:.2f} m/s, F = {F:.2f} N")

print("\nThe calculated coefficients for the equation are:")
print(f"  a = F = {a:.4e}")
print(f"  b = (0.5*M*v_0^2) - (F*l_0) - (G*m*M/l_0) = {b:.4e}")
print(f"  c = G*m*M = {c:.4e}")

print(f"\nSo, the final equation to solve is:")
print(f"({a:.4e}) * l_max^2 + ({b:.4e}) * l_max + ({c:.4e}) = 0")

# Calculate the discriminant to check for real solutions
discriminant = b**2 - 4 * a * c

if discriminant < 0:
    print("\nThe discriminant is negative. The spaceship's kinetic energy never reaches zero.")
    print("Therefore, the maximum distance is infinite.")
    l_max = float('inf')
else:
    # Use the quadratic formula to find the two roots
    sqrt_discriminant = math.sqrt(discriminant)
    root1 = (-b + sqrt_discriminant) / (2 * a)
    root2 = (-b - sqrt_discriminant) / (2 * a)

    print(f"\nThe two solutions for the distance are: {root1:.4e} m and {root2:.4e} m.")
    
    # The maximum distance is the first point (smallest distance > l_0) where the ship stops.
    if root1 > l_0 and root2 > l_0:
        l_max = min(root1, root2)
    elif root1 > l_0:
        l_max = root1
    elif root2 > l_0:
        l_max = root2
    else:
        l_max = float('nan') # Should not happen with these parameters

    print(f"\nThe physically relevant solution is the first point the spaceship reaches with zero velocity.")
    print(f"The maximum distance l_max is {l_max:.4e} meters.")

# Final answer in the requested format
# Convert the result to a float to remove scientific notation for the final answer block
final_answer = float(f"{l_max:.4e}")
<<<1384382.497554585>>>