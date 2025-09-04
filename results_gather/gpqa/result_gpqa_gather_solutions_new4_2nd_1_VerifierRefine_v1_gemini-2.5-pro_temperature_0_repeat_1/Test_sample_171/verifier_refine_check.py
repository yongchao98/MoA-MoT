import sympy
from sympy import symbols, Eq, exp, log, factor, together

def check_answer():
    """
    This function checks the correctness of the provided answer by symbolically deriving the relationship
    between the temperatures of the two stars based on the Boltzmann equation and the problem's conditions.
    """
    # 1. Define the symbolic variables
    # T1, T2 are temperatures, DeltaE is energy difference, k is Boltzmann constant.
    # C represents the constant ratio of statistical weights (g_j / g_i).
    # R1, R2 are the excitation ratios for star 1 and star 2.
    T1, T2, DeltaE, k, C, R1, R2 = symbols('T_1 T_2 DeltaE k C R_1 R_2', positive=True, real=True)

    # 2. State the Boltzmann Equation for each star
    # R = C * exp(-DeltaE / (k * T))
    eq_star1 = Eq(R1, C * exp(-DeltaE / (k * T1)))
    eq_star2 = Eq(R2, C * exp(-DeltaE / (k * T2)))

    # 3. Incorporate the problem's condition: "star_1 are twice as excited ... when compared to ... star_2"
    # This means R1 = 2 * R2
    condition = Eq(R1, 2 * R2)

    # 4. Substitute the Boltzmann expressions into the condition
    # This replaces R1 and R2 with their expressions in terms of temperature.
    subst_eq = condition.subs([(R1, eq_star1.rhs), (R2, eq_star2.rhs)])
    # -> C * exp(-DeltaE / (k * T1)) = 2 * C * exp(-DeltaE / (k * T2))

    # The constant C cancels out from both sides. Sympy handles this implicitly when solving.
    # We can simplify by dividing by C.
    simplified_eq = Eq(subst_eq.lhs / C, subst_eq.rhs / C)
    # -> exp(-DeltaE / (k * T1)) = 2 * exp(-DeltaE / (k * T2))

    # 5. Solve for the relationship between temperatures by taking the natural logarithm.
    # First, rearrange to isolate the constant '2'.
    rearranged_eq = Eq(simplified_eq.lhs / exp(-DeltaE / (k * T2)), 2)
    # -> exp(-DeltaE/(k*T1)) / exp(-DeltaE/(k*T2)) = 2
    # Using exponent rules, this is exp(DeltaE/(k*T2) - DeltaE/(k*T1)) = 2

    # Now, take the natural log of both sides.
    log_eq = Eq(log(rearranged_eq.lhs), log(rearranged_eq.rhs))
    # -> log(exp(...)) = log(2)
    # -> DeltaE/(k*T2) - DeltaE/(k*T1) = log(2)

    # 6. Perform algebraic manipulation to match the format of the options.
    # Factor out the common term (DeltaE / k).
    factored_eq_lhs = factor(log_eq.lhs)
    # -> (DeltaE/k) * (1/T2 - 1/T1)
    
    # Combine the fractions in the parenthesis.
    combined_fractions_lhs = together(factored_eq_lhs)
    # -> (DeltaE/k) * (T1 - T2) / (T1 * T2)

    # So, the full derived relationship is:
    # (DeltaE * (T1 - T2)) / (k * T1 * T2) = log(2)

    # 7. Use the given numerical values: DeltaE ≈ 1.38e-23 J and k ≈ 1.38e-23 J/K.
    # This means the ratio DeltaE / k is approximately 1.
    # We substitute this into our derived expression for the left-hand side.
    final_lhs = combined_fractions_lhs.subs(DeltaE / k, 1)
    # -> (T1 - T2) / (T1 * T2)

    # The derived equation is: log(2) = (T1 - T2) / (T1 * T2)

    # 8. Compare the derived result with the provided answer options.
    # The final answer given is 'D'. Let's check if our derivation matches option D.
    
    # Define the expressions for each option's right-hand side
    options = {
        'A': T2 / T1,
        'B': (T1 + T2) / (T1 * T2),
        'C': (T1 - T2) / (T1 * T2)**2,
        'D': (T1 - T2) / (T1 * T2)
    }
    
    chosen_answer_option = 'D'
    expression_from_answer = options[chosen_answer_option]

    # Check if the derived expression matches the expression from the chosen answer.
    if final_lhs == expression_from_answer:
        return "Correct"
    else:
        return (f"Incorrect. The derivation leads to the expression '{final_lhs}' for the right side of the equation. "
                f"The expression for the chosen answer '{chosen_answer_option}' is '{expression_from_answer}'. "
                f"These do not match.")

# Run the check
result = check_answer()
print(result)