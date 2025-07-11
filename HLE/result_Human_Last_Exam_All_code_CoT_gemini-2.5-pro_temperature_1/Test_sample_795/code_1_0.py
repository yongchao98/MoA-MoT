def print_magnetization_formula():
    """
    This function provides the analytical expression for the initial magnetization
    curve of a superconducting bar as described in the problem, using the
    Bean critical-state model for an infinite slab approximation.
    """

    # Explanation of the physical model and the result.
    intro_message = (
        "The problem is solved using the Bean critical-state model for a superconducting\n"
        "slab of width '2a' in a parallel applied magnetic field 'H'. The initial\n"
        "magnetization curve describes the behavior as H is increased from zero.\n\n"
        "The key result is that the magnetization 'M' has a quadratic dependence on 'H'\n"
        "until the field fully penetrates the superconductor. This occurs at a field\n"
        "known as the penetration field, H_p = Jc * a, where 'Jc' is the critical\n"
        "current density.\n"
    )

    # The analytical formula.
    formula = "M(H) = -H + H**2 / (2 * a * Jc)"

    # The range of validity for the formula.
    validity = "This expression is valid for the range of applied fields: 0 <= H <= a * Jc."

    print(intro_message)
    print("The analytical expression for the initial magnetization is:")
    print(formula)
    print(validity)
    
    # Per the instructions, we now explicitly output the numbers in the final equation.
    # The general form of the equation is: M = (c1 * H) + (c2 * H**2) / (c3 * a * Jc)
    c1 = -1
    c2 = 1
    c3 = 2
    
    print("\n" + "="*50 + "\n")
    print("To meet the output requirements, here are the numerical coefficients from the equation:")
    print("M = (" + str(c1) + " * H) + (" + str(c2) + " * H**2) / (" + str(c3) + " * a * Jc)")
    print(f"\n- The numerical coefficient for the linear term (H) is: {c1}")
    print(f"- The numerical coefficient for the quadratic term (H**2) in the numerator is: {c2}")
    print(f"- The numerical coefficient in the denominator is: {c3}")


if __name__ == "__main__":
    print_magnetization_formula()
