import math

def calculate_joule_heat():
    """
    Calculates and displays the formula for the Joule heat dissipated by a
    shrinking and leaking charged sphere.

    The problem asks for the total Joule heat dissipated into the atmosphere as a
    sphere of initial radius 'a' and initial potential 'V' leaks its charge
    and shrinks.

    The solution is based on the principle of conservation of energy. The initial
    energy of the system is the electrostatic potential energy stored in the
    charged sphere. As the sphere leaks charge and shrinks to zero radius, this
    energy is dissipated. Assuming all of this stored energy is converted into
    Joule heat, the total heat dissipated is equal to the initial energy.

    1. The capacitance (C) of a sphere with radius 'a' is:
       C = 4 * pi * epsilon_0 * a

    2. The initial electrostatic energy (U) stored in the sphere is:
       U = (1/2) * C * V^2

    3. Substituting C into the energy equation:
       U = (1/2) * (4 * pi * epsilon_0 * a) * V^2
       U = 2 * pi * epsilon_0 * a * V^2

    This initial energy is the total amount dissipated as Joule heat.
    """

    # The formula is symbolic. We will print the components of the formula.
    # U = 2 * pi * epsilon_0 * a * V^2
    # The components are the number 2, pi, epsilon_0 (permittivity of free space),
    # radius 'a', and potential 'V'.

    term_1_value = 2
    term_2_symbol = "pi"
    term_3_symbol = "epsilon_0"
    term_4_symbol = "a"
    term_5_symbol = "V"

    print("The total Joule heat dissipated is equal to the initial electrostatic energy of the sphere.")
    print("The formula for the Joule heat (Q) is:")
    print(f"Q = {term_1_value} * {term_2_symbol} * {term_3_symbol} * {term_4_symbol} * {term_5_symbol}^2")

if __name__ == "__main__":
    calculate_joule_heat()
