import math

# --- Define Constants and Variables ---
# I will use a set of parameters that ensures a finite maximum distance exists.
# This typically requires the initial gravitational force to be greater than the
# applied force, and the initial velocity to be within a specific range.

# Gravitational constant in N m^2 / kg^2
G = 6.674e-11
# Mass of the asteroid A in kg
m = 1.0e15
# Mass of the spaceship B in kg
M = 5.0e4
# Initial distance between A and B in meters
l0 = 8000.0
# Initial speed of the spaceship in m/s
v0 = 0.5
# Additional force applied to the spaceship in Newtons
F = 30.0

# --- Physics Calculation using the Work-Energy Theorem ---

# The work-energy theorem leads to the equation:
# F*l_max + (G*M*m)/l_max = F*l0 + (G*M*m)/l0 - 0.5*M*v0^2
#
# This can be rearranged into a quadratic equation for l_max:
# (F) * l_max^2 - (F*l0 + G*M*m/l0 - 0.5*M*v0^2) * l_max + (G*M*m) = 0
# This is in the form a*x^2 + b*x + c = 0, where x = l_max.

# Calculate the terms
G_M_m = G * M * m
initial_kinetic_energy = 0.5 * M * v0**2
energy_constant_K = F * l0 + G_M_m / l0 - initial_kinetic_energy

# --- Solve the Quadratic Equation ---

# Coefficients of the quadratic equation a*x^2 + b*x + c = 0
a = F
b = -energy_constant_K
c = G_M_m

# Calculate the discriminant (b^2 - 4ac)
discriminant = b**2 - 4 * a * c

l_max = 0
# Ensure a real solution exists
if discriminant >= 0:
    # Calculate the two roots
    sqrt_discriminant = math.sqrt(discriminant)
    root1 = (-b + sqrt_discriminant) / (2 * a)
    root2 = (-b - sqrt_discriminant) / (2 * a)

    # The spaceship starts at l0 and moves outwards. The maximum distance is the
    # first turning point it reaches, which is the smaller of the two roots
    # that is greater than l0.
    if root1 > l0 and root2 > l0:
        l_max = min(root1, root2)
    elif root1 > l0:
        l_max = root1
    elif root2 > l0:
        l_max = root2
    else:
        # This case implies no turning point is reached beyond l0
        l_max = float('inf') 
else:
    # A negative discriminant means the ship's energy is not sufficient to
    # reach a turning point; it escapes to infinity.
    l_max = float('inf')

# --- Output the Result ---

print("Based on the work-energy theorem, we solve the following quadratic equation for l_max:")
# Using integer values for display clarity where appropriate
print(f"({F:.1f}) * l_max^2 - ({energy_constant_K:,.2f}) * l_max + ({G_M_m:,.2f}) = 0\n")

if l_max != float('inf'):
    print(f"The calculated maximum distance l_max is: {l_max:,.2f} meters.")
else:
    print("The spaceship has enough energy to escape to infinity; no finite maximum distance is reached.")

# For the final answer format
final_answer = l_max

# Print the final numerical answer for extraction
# print(f"\n<<<__{final_answer:.2f}__>>>")
# The format requested is <<<answer content>>> without extra text. Let's make sure it's correct.
print(f"\n<<<{final_answer:.2f}>>>")
