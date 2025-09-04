import sympy as sp

def check_stellar_temperature_equation():
    """
    Verifies the correct relationship between stellar temperatures T1 and T2
    based on the problem description using symbolic mathematics.
    """
    # Define symbolic variables for the physical quantities.
    # T1 and T2 are positive temperatures.
    T1, T2, delta_E, k = sp.symbols('T_1 T_2 delta_E k', positive=True)

    # The problem states that the excitation in star 1 is twice that in star 2.
    # Based on the Boltzmann equation, R(T) = C * exp(-delta_E / (k*T)),
    # the relationship R(T1) = 2 * R(T2) simplifies to:
    # exp(-delta_E / (k*T1)) = 2 * exp(-delta_E / (k*T2))
    
    # We create this equation symbolically.
    initial_equation = sp.Eq(sp.exp(-delta_E / (k * T1)), 2 * sp.exp(-delta_E / (k * T2)))

    # To solve for the relationship, we take the natural logarithm of both sides.
    # Using sympy.log and the .expand() method correctly applies log rules.
    log_equation = sp.Eq(sp.log(initial_equation.lhs), sp.log(initial_equation.rhs)).expand(force=True)
    
    # The resulting equation is: -delta_E/(k*T1) = log(2) - delta_E/(k*T2)
    # Now, we solve this equation for log(2) to get the expression we need to check.
    derived_expression_for_ln2 = sp.solve(log_equation, sp.log(2))[0]

    # The result is delta_E/(k*T2) - delta_E/(k*T1).
    # We simplify this by finding a common denominator.
    simplified_derived_expr = sp.simplify(derived_expression_for_ln2)

    # The problem is designed such that the ratio delta_E / k is approximately 1.
    # We substitute this simplification into our derived formula.
    final_derived_expr = simplified_derived_expr.subs(delta_E / k, 1)

    # The final derived formula is: (T1 - T2) / (T1 * T2)

    # Now, we define the expressions from the multiple-choice options.
    # The final answer to check is 'A'.
    options = {
        'A': (T1 - T2) / (T1 * T2),
        'B': T2 / T1,
        'C': (T1 + T2) / (T1 * T2),
        'D': (T1 - T2) / (T1 * T2)**2
    }
    
    answer_to_check = 'A'
    expression_from_answer = options[answer_to_check]

    # Finally, we check if our rigorously derived expression matches the expression
    # from the provided answer 'A'.
    if final_derived_expr == expression_from_answer:
        return "Correct"
    else:
        return (f"Incorrect. The derivation from first principles yields ln(2) = {final_derived_expr}. "
                f"The provided answer '{answer_to_check}' corresponds to the expression {expression_from_answer}, "
                f"which does not match the derivation.")

# Execute the check
result = check_stellar_temperature_equation()
print(result)