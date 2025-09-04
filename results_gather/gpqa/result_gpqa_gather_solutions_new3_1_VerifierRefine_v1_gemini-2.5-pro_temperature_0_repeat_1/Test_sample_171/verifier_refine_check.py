import sympy

def check_stellar_temperature_equation():
    """
    Symbolically derives the relationship between the temperatures of the two stars
    based on the problem description and checks it against the provided answer.
    """
    # 1. Define all the symbols involved in the problem.
    # T_1, T_2: Temperatures of the two stars
    # delta_E: Energy difference between the atomic levels
    # k: Boltzmann constant
    # g_ratio: The ratio of statistical weights (g_j / g_i), which is a constant
    T_1, T_2, delta_E, k, g_ratio = sympy.symbols('T_1 T_2 delta_E k g_ratio', positive=True)

    # 2. Write the Boltzmann equation for the excitation ratio (R) for each star.
    # R = g_ratio * exp(-delta_E / (k * T))
    R_1 = g_ratio * sympy.exp(-delta_E / (k * T_1))
    R_2 = g_ratio * sympy.exp(-delta_E / (k * T_2))

    # 3. Apply the core constraint from the question: "star_1 is twice as excited as star_2".
    # This means R_1 = 2 * R_2.
    # We create an equation object representing this relationship.
    main_equation = sympy.Eq(R_1, 2 * R_2)

    # 4. Simplify the equation. The 'g_ratio' term appears on both sides and can be canceled
    # as long as it's not zero, which is physically true.
    # sympy.cancel() can do this, or we can divide both sides by g_ratio.
    simplified_equation = sympy.Eq(main_equation.lhs / g_ratio, main_equation.rhs / g_ratio)

    # 5. To solve for the temperatures, we take the natural logarithm (ln) of both sides.
    # This helps to remove the exponentials.
    log_equation = sympy.Eq(sympy.log(simplified_equation.lhs), sympy.log(simplified_equation.rhs))

    # 6. Use logarithm properties to expand the terms. ln(e^x) = x and ln(a*b) = ln(a) + ln(b).
    # sympy.expand_log() or sympy.logexpand() can do this.
    expanded_log_equation = sympy.expand_log(log_equation, force=True)

    # 7. Rearrange the equation to solve for ln(2).
    # The equation is currently: -delta_E/(k*T_1) = ln(2) - delta_E/(k*T_2)
    # We want to isolate ln(2).
    try:
        derived_rhs = sympy.solve(expanded_log_equation, sympy.log(2))[0]
    except IndexError:
        return "Failed to solve for ln(2) during derivation."

    # 8. The problem is designed such that the ratio delta_E / k is approximately 1.
    # We substitute this key simplification into our derived expression.
    # The derived expression is currently (delta_E/k) * (1/T_2 - 1/T_1).
    final_derived_expr = derived_rhs.subs(delta_E / k, 1)

    # 9. Now, we define the expression from the proposed correct answer, Option C.
    # Option C is: ln(2) = [ (T_1 - T_2) / (T1*T2)]
    option_c_expr = (T_1 - T_2) / (T_1 * T_2)

    # 10. Finally, check if our derived expression is mathematically equivalent to Option C's expression.
    # sympy.simplify() can help ensure both are in a standard form for comparison.
    if sympy.simplify(final_derived_expr - option_c_expr) == 0:
        return "Correct"
    else:
        return (f"Incorrect. The derivation leads to the expression: {final_derived_expr}. "
                f"The expression for option C is: {option_c_expr}. These are not equivalent.")

# Run the check
result = check_stellar_temperature_equation()
print(result)