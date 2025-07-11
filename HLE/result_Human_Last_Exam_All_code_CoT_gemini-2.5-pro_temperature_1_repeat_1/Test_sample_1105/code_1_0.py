import math

def calculate_max_distance():
    """
    Calculates the maximum distance a spaceship can achieve from an asteroid.

    The calculation is based on the work-energy theorem, which leads to a
    quadratic equation for the maximum distance, l_max.
    """
    # Define the physical constants and variables in SI units.
    # Mass of the asteroid A
    m = 1.0e18  # kg
    # Mass of the spaceship B
    M = 5.0e5   # kg
    # Initial distance between A and B
    l_0 = 2.0e5 # m
    # Initial speed of the spaceship
    v_0 = 150.0 # m/s
    # Additional constant force applied to the spaceship
    F = 2000.0  # N
    # Gravitational constant
    G = 6.674e-11 # N m^2 / kg^2

    # The governing equation is derived from the work-energy principle:
    # Work_by_F + Work_by_Gravity = Change_in_Kinetic_Energy
    # F*(l_max - l_0) + G*m*M*(1/l_max - 1/l_0) = 0 - 1/2*M*v_0^2
    # This can be rearranged to an energy conservation form:
    # F*l_max + G*m*M/l_max = F*l_0 + G*m*M/l_0 - 1/2*M*v_0^2

    print("The governing equation based on the work-energy theorem is:")
    print("F * l_max + (G * m * M) / l_max = F * l_0 + (G * m * M) / l_0 - 0.5 * M * v_0^2\n")

    print("Substituting the given values into the equation:")
    print(f"{F} * l_max + ({G} * {m} * {M}) / l_max = {F} * {l_0} + ({G} * {m} * {M}) / {l_0} - 0.5 * {M} * {v_0}^2\n")

    # To solve for l_max, we form a quadratic equation: a*x^2 + b*x + c = 0, where x = l_max.
    # F * l_max^2 - (F*l_0 + G*m*M/l_0 - 0.5*M*v_0^2) * l_max + G*m*M = 0

    # Calculate the right-hand side of the energy equation, which is a constant.
    rhs = F * l_0 + (G * m * M) / l_0 - 0.5 * M * v_0**2

    # Define the coefficients of the quadratic equation.
    a = F
    b = -rhs
    c = G * m * M

    # Calculate the discriminant to find the roots.
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        result_message = "No real solution exists. The spaceship cannot reach a point of zero velocity under these conditions."
        final_answer = "Error: No real solution"
    else:
        # The larger root corresponds to the maximum distance reached after moving away from l_0.
        l_max = (-b + math.sqrt(discriminant)) / (2 * a)
        result_message = f"The maximum distance l_max is: {l_max:.4e} meters."
        final_answer = l_max

    print(result_message)
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    calculate_max_distance()