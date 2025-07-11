import math

# --- Step 1: Define physical constants and initial conditions ---
# These values are chosen to be physically consistent for a solvable scenario.

# Gravitational constant in N * m^2 / kg^2
G = 6.67430e-11 
# Mass of the asteroid A in kg
m = 1.0e15
# Mass of the spaceship B in kg
M = 5.0e4
# Initial distance between A and B in meters
l_0 = 5000.0
# Initial speed of the spaceship B in m/s
v_0 = 8.0
# Additional constant force applied to the spaceship B in Newtons
F = 110.0

# --- Step 2: Formulate the quadratic equation for l_max ---
# The work-energy theorem leads to the equation:
# F*l_max + (G*M*m)/l_max = F*l_0 + (G*M*m)/l_0 - 0.5*M*v_0^2
#
# Let C = F*l_0 + (G*M*m)/l_0 - 0.5*M*v_0^2
# The equation becomes: F*l_max + (G*M*m)/l_max = C
#
# Rearranging into the quadratic form a*x^2 + b*x + c = 0, where x = l_max:
# F*(l_max)^2 - C*l_max + (G*M*m) = 0

# Calculate the constant term C based on initial conditions
C = (F * l_0) + (G * M * m / l_0) - (0.5 * M * v_0**2)

# Define the coefficients of the quadratic equation a*x^2 + b*x + c = 0
a = F
b = -C
c = G * M * m

# --- Step 3: Solve the quadratic equation ---

# Calculate the discriminant (delta = b^2 - 4ac)
delta = b**2 - 4 * a * c

# Check for a valid physical solution
if delta < 0:
    print("No real solution exists. The spaceship either escapes to infinity or cannot move away under these conditions.")
else:
    # Calculate the two roots of the equation
    sqrt_delta = math.sqrt(delta)
    root1 = (-b + sqrt_delta) / (2 * a)
    root2 = (-b - sqrt_delta) / (2 * a)
    
    # The maximum distance l_max is the larger of the two roots.
    # The smaller root corresponds to the initial position l_0.
    l_max = max(root1, root2)
    
    # --- Step 4: Output the results ---
    
    # Print the equation with its numerical coefficients
    print("The final equation for the maximum distance (l_max) is a quadratic equation:")
    print(f"({a:.2f}) * l_max^2 + ({b:.2f}) * l_max + ({c:.2e}) = 0")
    print("\nSolving this equation gives the maximum distance.")
    
    # Print the final calculated maximum distance
    print(f"\nThe initial distance was: {l_0:.2f} m")
    print(f"The maximum distance the spaceship can achieve is: {l_max:.2f} m")
    
    print("\n<<<" + str(l_max) + ">>>")
