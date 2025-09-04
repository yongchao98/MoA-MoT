import sympy

def check_answer():
    """
    Checks the correctness of the provided Liénard-Wiechert potentials.

    The function performs two checks:
    1. Direct Comparison: It compares the expressions from answer D with the standard,
       textbook formulas for the Liénard-Wiechert potentials.
    2. Consistency Check: It verifies the known relationship between the scalar and
       vector potentials, A = (v/c^2) * V, using the expressions from answer D and
       the fundamental relation c^2 = 1/(epsilon_0 * mu_0).
    """
    # Define the symbols used in the equations
    # Scalars
    q = sympy.Symbol('q')          # charge
    c = sympy.Symbol('c', positive=True)          # speed of light
    epsilon_0 = sympy.Symbol('epsilon_0') # permittivity of free space
    mu_0 = sympy.Symbol('mu_0')       # permeability of free space
    d = sympy.Symbol('d', positive=True)          # magnitude of vector d
    pi = sympy.pi

    # To handle vector quantities symbolically, we represent the dot product
    # d_vec . v_vec as a single scalar symbol.
    d_dot_v = sympy.Symbol('d_dot_v')

    # The vector v_vec appears linearly, so we can treat it as a symbolic constant.
    v_vec = sympy.Symbol('v_vec')

    # --- Expressions from the chosen answer (D) ---
    # Denominator term, common to both V and A in option D
    denominator_D = 4 * pi * (d * c - d_dot_v)

    # Scalar Potential V from option D
    V_D = (q * c) / (epsilon_0 * denominator_D)

    # Vector Potential A from option D
    A_D = (mu_0 * q * c * v_vec) / denominator_D

    # --- Standard (Correct) Liénard-Wiechert Potentials ---
    # These are the known correct formulas from physics textbooks.
    V_correct = (q * c) / (4 * pi * epsilon_0 * (d * c - d_dot_v))
    A_correct = (mu_0 * q * c * v_vec) / (4 * pi * (d * c - d_dot_v))

    # --- Check 1: Direct comparison with standard formulas ---
    # We use sympy.simplify to see if the difference between the answer's expression
    # and the correct expression is zero.
    if sympy.simplify(V_D - V_correct) != 0:
        return f"Incorrect: The scalar potential V is wrong. Expected {V_correct} but got {V_D}."

    if sympy.simplify(A_D - A_correct) != 0:
        return f"Incorrect: The vector potential A is wrong. Expected {A_correct} but got {A_D}."

    # --- Check 2: Internal consistency check ---
    # The Liénard-Wiechert potentials must satisfy the relation A = (v/c^2) * V.
    # Let's verify if the expressions from option D satisfy this.
    # We also use the fundamental relation c^2 = 1 / (epsilon_0 * mu_0).

    # Calculate (v/c^2) * V using the expression for V from option D
    V_times_v_div_c2 = (v_vec / c**2) * V_D

    # Now, let's see if this equals the expression for A from option D.
    # We substitute mu_0 = 1 / (epsilon_0 * c**2) into the expression for A_D.
    A_D_substituted = A_D.subs(mu_0, 1 / (epsilon_0 * c**2))

    # The expressions should be identical.
    if sympy.simplify(A_D_substituted - V_times_v_div_c2) != 0:
        return (f"Incorrect: The potentials in answer D are not self-consistent. "
                f"The relation A = (v/c^2) * V is not satisfied. "
                f"Based on V, A should be {V_times_v_div_c2}, but the answer gives an expression "
                f"equivalent to {A_D_substituted}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)