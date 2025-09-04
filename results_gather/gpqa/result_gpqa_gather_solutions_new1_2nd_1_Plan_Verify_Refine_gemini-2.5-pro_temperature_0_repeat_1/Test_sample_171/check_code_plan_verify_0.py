import sympy
from sympy import symbols, Eq, exp, log, simplify

def check_stellar_temperature_equation():
    """
    This function verifies the correctness of the answer by symbolically deriving
    the relationship between the temperatures of the two stars based on the
    problem description.

    The derivation follows these steps:
    1.  Define the Boltzmann equation for the excitation ratio for each star.
    2.  Apply the condition that the excitation ratio of star_1 is twice that of star_2.
    3.  Solve the resulting equation for ln(2).
    4.  Apply the simplification that the ratio of the energy difference (ΔE) to the
        Boltzmann constant (k) is approximately 1.
    5.  Compare the derived expression to the expression from the chosen answer.
    """
    try:
        # 1. Define symbolic variables based on the problem statement.
        # T1, T2: temperatures of star_1 and star_2
        # delta_E: energy difference between the energy levels
        # k: Boltzmann constant
        # C: ratio of statistical weights (g_j / g_i), a constant for the transition
        T1, T2, delta_E, k, C = symbols('T1 T2 delta_E k C', positive=True, real=True)

        # 2. State the governing principle: The Boltzmann Equation for the excitation ratio R.
        # R(T) = C * exp(-delta_E / (k * T))
        R1 = C * exp(-delta_E / (k * T1))
        R2 = C * exp(-delta_E / (k * T2))

        # 3. Apply the problem's central condition: R1 = 2 * R2
        condition_eq = Eq(R1, 2 * R2)

        # 4. Solve the equation for ln(2).
        # First, take the natural logarithm of both sides.
        # sympy.log is the natural logarithm (ln).
        log_eq = Eq(log(condition_eq.lhs), log(condition_eq.rhs))

        # Expand the logarithms using log rules: log(a*b) = log(a) + log(b)
        # and log(exp(x)) = x.
        # log(C) - delta_E / (k * T1) = log(2) + log(C) - delta_E / (k * T2)
        # Rearranging for log(2) gives:
        # log(2) = delta_E / (k * T2) - delta_E / (k * T1)
        
        # Let's have sympy solve for log(2) to be robust.
        # We can isolate log(2) by rearranging the expanded log equation.
        derived_log2_expr = sympy.solve(log_eq, log(2))
        
        # The solver returns a list, so we take the first element.
        if not derived_log2_expr:
            return "Error: Could not solve for log(2) symbolically."
        derived_expr_with_constants = derived_log2_expr[0]

        # 5. Apply the crucial simplification from the problem's numerical values:
        # delta_E ≈ 1.38e-23 J and k ≈ 1.38e-23 J/K, so delta_E / k ≈ 1.
        # We substitute the ratio (delta_E / k) with 1.
        final_derived_expr = derived_expr_with_constants.subs(delta_E / k, 1)
        
        # Simplify the resulting expression.
        final_derived_expr = simplify(final_derived_expr)

        # 6. Define the expression from the chosen answer (B).
        # The final answer is <<<B>>>.
        # The options listed in the final answer text are:
        # A) ln(2) = [ (T_1 + T_2) / (T1*T2)]
        # B) ln(2) = [ (T_1 - T_2) / (T1*T2)]
        # C) ln(2) = [ (T_1 - T_2) / (T1*T2)^2 ]
        # D) ln(2) = [ T_2 / T_1]
        # Answer B corresponds to the expression (T1 - T2) / (T1 * T2).
        answer_b_expr = (T1 - T2) / (T1 * T2)

        # 7. Compare the derived expression with the expression from the answer.
        # The .equals() method is a robust way to check for symbolic equality.
        if final_derived_expr.equals(answer_b_expr):
            return "Correct"
        else:
            return (f"Incorrect. The derivation from first principles results in the expression for ln(2) being: "
                    f"{final_derived_expr}. The expression from the chosen answer B is: {answer_b_expr}. "
                    f"These two expressions are not equivalent.")

    except Exception as e:
        return f"An error occurred during the checking process: {e}"

# Run the check and print the result.
print(check_stellar_temperature_equation())