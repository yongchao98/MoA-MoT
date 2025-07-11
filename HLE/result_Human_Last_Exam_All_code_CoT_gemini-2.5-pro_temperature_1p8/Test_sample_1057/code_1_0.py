def solve_sphere_joule_heat():
    """
    This script presents the derived formula for the Joule heat dissipated
    by a shrinking and charge-leaking sphere.
    """
    # The final formula for Joule heat (H) is of the form:
    # H = (numerator/denominator) * pi * epsilon_0 * a * V^2
    # Based on the physics derivation:
    # Initial stored energy: U_initial = 2 * pi * epsilon_0 * a * V^2
    # Mechanical work done on the sphere against electrostatic repulsion:
    # W_mech = (2/3) * pi * epsilon_0 * a * V^2
    # Total heat H = U_initial + W_mech

    # Numerator of the final coefficient (2 + 2/3 = 8/3)
    numerator = 8
    # Denominator of the final coefficient
    denominator = 3

    print("The total Joule heat (H) dissipated into the atmosphere is given by the equation:")
    # Using Unicode for mathematical symbols: π (pi), ε₀ (epsilon_naught), ² (superscript two)
    print(f"H = ({numerator}/{denominator}) * \u03C0 * \u03B5\u2080 * a * V\u00B2")

solve_sphere_joule_heat()