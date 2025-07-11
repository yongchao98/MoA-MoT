import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid
    based on initial conditions and applied forces.
    """
    # Physical constants and problem parameters
    # Gravitational constant in N(m/kg)^2
    G = 6.67430e-11
    # Mass of the asteroid (m) in kg
    m = 6e20
    # Mass of the spaceship (M) in kg
    M = 5e4
    # Initial distance (l_0) in meters
    l_0 = 10000.0
    # Initial speed (v_0) in m/s
    v_0 = 100.0
    # Additional applied force (F) in Newtons
    F = 1e7

    # Print the given parameters
    print("Given Parameters:")
    print(f"  Gravitational Constant (G): {G} N(m/kg)^2")
    print(f"  Mass of Asteroid (m): {m:.1e} kg")
    print(f"  Mass of Spaceship (M): {M:.1e} kg")
    print(f"  Initial Distance (l_0): {l_0} m")
    print(f"  Initial Speed (v_0): {v_0} m/s")
    print(f"  Applied Force (F): {F:.1e} N\n")

    # The problem is solved by setting up an energy balance equation, which
    # can be rearranged into a quadratic equation for l_max:
    # a*l_max^2 + b*l_max + c = 0

    # Calculate the coefficients of the quadratic equation
    # a = F
    # b = -(F*l_0 + G*m*M/l_0 - 0.5*M*v_0**2)
    # c = G*m*M

    a_coeff = F
    b_coeff = -(F * l_0 + (G * m * M) / l_0 - 0.5 * M * v_0**2)
    c_coeff = G * m * M

    # Print the quadratic equation with its numerical coefficients
    print("The problem reduces to solving the quadratic equation for l_max:")
    print(f"({a_coeff:.4e}) * l_max^2 + ({b_coeff:.4e}) * l_max + ({c_coeff:.4e}) = 0\n")

    # Calculate the discriminant
    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff

    if discriminant < 0:
        print("No real solution exists. The initial energy and applied force are not sufficient for the spaceship to move further away.")
    else:
        # Calculate the two roots using the quadratic formula
        sqrt_discriminant = math.sqrt(discriminant)
        l_max1 = (-b_coeff + sqrt_discriminant) / (2 * a_coeff)
        l_max2 = (-b_coeff - sqrt_discriminant) / (2 * a_coeff)

        # The maximum distance is the larger of the two roots, as the spaceship
        # starts at l_0 and moves outwards to its farthest point.
        l_max = max(l_max1, l_max2)

        print(f"The calculated maximum distance (l_max) is: {l_max:.4f} meters")
        # Return the value for the final answer format
        return l_max

# Run the calculation and store the result
final_answer = calculate_max_distance()

# The final answer in the required format
# To avoid printing the answer twice, we only output the formatted string.
# We use a placeholder check to ensure the function ran correctly.
if final_answer is not None:
    print(f"<<<{final_answer:.4f}>>>")
