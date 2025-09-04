import sympy

def check_lienard_wiechert_potentials():
    """
    Symbolically verifies the expressions for the Liénard-Wiechert potentials.

    This function checks if the expressions in Option B match the standard
    textbook formulas for the scalar and vector potentials of a moving point charge.
    It also verifies the relationship A = (v/c^2)V.
    """
    # Define the symbols used in the equations.
    # We assume they are real and positive where appropriate.
    q, c, epsilon_0, mu_0, d, pi = sympy.symbols('q c epsilon_0 mu_0 d pi', real=True, positive=True)
    
    # For vector quantities, we use symbols to represent the final scalar or vector results.
    # v_vec represents the vector v.
    # d_dot_v represents the scalar result of the dot product (d . v).
    v_vec = sympy.Symbol('v_vec') 
    d_dot_v = sympy.Symbol('d_dot_v')

    # --- Part 1: Define the expressions from the given answer (Option B) ---
    V_answer = (q * c) / (4 * pi * epsilon_0 * (d * c - d_dot_v))
    A_answer = (mu_0 * q * c * v_vec) / (4 * pi * (d * c - d_dot_v))

    # --- Part 2: Define the standard Liénard-Wiechert potential formulas ---
    # The standard denominator is R(1 - (R_hat . v)/c), which simplifies to (dc - d.v)/c
    # in the problem's notation.
    denominator_factor = (d * c - d_dot_v) / c
    
    V_standard = (1 / (4 * pi * epsilon_0)) * (q / denominator_factor)
    A_standard = (mu_0 / (4 * pi)) * (q * v_vec / denominator_factor)

    # --- Part 3: Perform the checks ---

    # Check 1: Compare the scalar potential V from the answer with the standard formula.
    # We simplify the difference. If it's zero, they are equivalent.
    if sympy.simplify(V_answer - V_standard) != 0:
        return (f"Incorrect. The scalar potential V is wrong. "
                f"The expression from the answer simplifies to {sympy.simplify(V_answer)}, "
                f"while the standard expression simplifies to {sympy.simplify(V_standard)}.")

    # Check 2: Compare the vector potential A from the answer with the standard formula.
    if sympy.simplify(A_answer - A_standard) != 0:
        return (f"Incorrect. The vector potential A is wrong. "
                f"The expression from the answer simplifies to {sympy.simplify(A_answer)}, "
                f"while the standard expression simplifies to {sympy.simplify(A_standard)}.")

    # Check 3: Verify the relationship A = (v/c^2)V using the answer's expressions.
    # This is a known identity for Liénard-Wiechert potentials.
    A_derived_from_V = (v_vec / c**2) * V_answer
    
    # To compare A_derived_from_V with A_answer, we use the relation c^2 = 1/(epsilon_0 * mu_0).
    # We substitute mu_0 in A_answer to express it in terms of c and epsilon_0.
    difference = A_answer.subs(mu_0, 1 / (c**2 * epsilon_0)) - A_derived_from_V
    if sympy.simplify(difference) != 0:
        return (f"Incorrect. The relationship A = (v/c^2)V is not satisfied. "
                f"The difference A - (v/c^2)V simplifies to {sympy.simplify(difference)}, not zero.")

    # Check 4: Briefly analyze other options.
    # Option C has the wrong sign in the denominator.
    V_C = (q * c) / (4 * pi * epsilon_0 * (d * c + d_dot_v))
    if sympy.simplify(V_C - V_standard) == 0:
        # This should not happen, but it's a sanity check for our logic.
        return "Checker error: Option C with '+' sign evaluated as correct."
    
    # Options A and D are quasi-static approximations, not the general retarded potentials.
    # Their form is fundamentally different (e.g., denominator is just 'r' or 'd').
    # They are incorrect because they don't account for the (dc - d.v) retardation factor.

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_lienard_wiechert_potentials()
print(result)