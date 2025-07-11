import math

# Step 1: Define physical constants and initial parameters.
# Using values that lead to a finite maximum distance for demonstration.
# Gravitational constant in N * m^2 / kg^2
G = 6.67430e-11
# Mass of the asteroid A in kg
m_asteroid = 5.972e22  # A very large asteroid, or a small moon
# Mass of the spaceship B in kg
M_spaceship = 50000.0
# Initial distance between A and B in meters
l_0 = 1000000.0  # 1000 km
# Initial speed of the spaceship in m/s
v_0 = 200.0
# Additional applied force in Newtons
F = 1000.0

# Step 2: Formulate the quadratic equation for l_max:
# a*l_max^2 + b*l_max + c = 0
# derived from the work-energy theorem.

# The equation is: F * l_max^2 - C * l_max + (G*m*M) = 0
# where C = F*l_0 + (G*m*M)/l_0 - 0.5*M*v_0^2

# Calculate the components of the equation
GmM = G * m_asteroid * M_spaceship
initial_K = 0.5 * M_spaceship * v_0**2
constant_C = F * l_0 + GmM / l_0 - initial_K

# Coefficients of the quadratic equation a*x^2 + b*x + c = 0 where x = l_max
a = F
b = -constant_C
c = GmM

# Step 3: Solve the quadratic equation.
# First, calculate the discriminant (D = b^2 - 4ac)
discriminant = b**2 - 4 * a * c

print("The maximum distance, l_max, is a solution to the quadratic equation derived from the work-energy theorem.")
print("The final equation form is: a * (l_max)^2 + b * l_max + c = 0\n")
print("Based on the provided parameters, the equation with numerical values is:")
print(f"{a:.1f} * (l_max)^2 + ({b:.1f}) * l_max + ({c:.1e}) = 0\n")


# Check if a real solution exists (i.e., if the spaceship ever stops)
if discriminant < 0:
    print("The discriminant is negative. This means the spaceship has enough energy to overcome the gravitational pull and will never stop.")
    print("It will travel to an infinite distance, so a finite l_max does not exist.")
else:
    # Use the quadratic formula to find the two roots for the distance.
    # We choose the '+' sign to get the larger root, which represents the maximum distance reached.
    l_max = (-b + math.sqrt(discriminant)) / (2 * a)

    print("Solving for l_max (the larger root of the equation) gives the maximum distance.")
    print(f"The maximum distance achieved by the spaceship is: {l_max:,.2f} meters.")
    print(f"This is approximately {l_max/1000:,.2f} kilometers.")

    # The final answer in the requested format
    print("\n<<<", end="")
    print(f"{a} * l_max^2 - {constant_C} * l_max + {GmM} = 0", end="")
    print(">>>")
