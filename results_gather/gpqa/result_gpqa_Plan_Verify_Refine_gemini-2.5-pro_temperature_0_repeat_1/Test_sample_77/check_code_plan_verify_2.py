import sympy

def check_lienard_wiechert_potentials():
    """
    Checks the correctness of the given Liénard-Wiechert potentials.
    """
    # Define the symbols used in the physics problem.
    # We treat them as positive real numbers where applicable.
    q, c, epsilon_0, mu_0, d = sympy.symbols('q c epsilon_0 mu_0 d', real=True, positive=True)
    
    # The dot product d_vec . v_vec is a scalar.
    d_dot_v = sympy.Symbol('d_dot_v', real=True)
    
    # The velocity v_vec is a vector. We represent it as a symbol.
    # In these equations, it acts as a vector scaling factor.
    v_vec = sympy.Symbol('v_vec')

    # --- Part 1: Define the correct, standard Liénard-Wiechert potentials ---
    # These are the accepted formulas from electromagnetism textbooks.
    # The common denominator is (d*c - d_vec . v_vec)
    denominator_term = (d * c - d_dot_v)
    
    # Standard formula for scalar potential V
    V_correct = (q * c) / (4 * sympy.pi * epsilon_0 * denominator_term)
    
    # Standard formula for vector potential A
    A_correct = (mu_0 * q * c * v_vec) / (4 * sympy.pi * denominator_term)

    # --- Part 2: Define the potentials from the proposed answer (Option D) ---
    V_D = (q * c) / (4 * sympy.pi * epsilon_0 * (d * c - d_dot_v))
    A_D = (mu_0 * q * c * v_vec) / (4 * sympy.pi * (d * c - d_dot_v))

    # --- Part 3: Perform the verification ---

    # Check 1: Is the scalar potential V from Option D correct?
    # We simplify the difference. If it's zero, they are algebraically equivalent.
    if sympy.simplify(V_D - V_correct) != 0:
        return (f"Incorrect: The scalar potential V in option D is wrong. "
                f"The expression simplifies to {V_D}, but the correct expression is {V_correct}.")

    # Check 2: Is the vector potential A from Option D correct?
    if sympy.simplify(A_D - A_correct) != 0:
        return (f"Incorrect: The vector potential A in option D is wrong. "
                f"The expression simplifies to {A_D}, but the correct expression is {A_correct}.")

    # Check 3: Does the relationship A = (v/c^2)V hold for Option D?
    # This is a fundamental property of Liénard-Wiechert potentials,
    # derived from the relation c^2 = 1 / (mu_0 * epsilon_0).
    # A/V should equal mu_0 * epsilon_0 * v_vec
    
    # Calculate the ratio from Option D's formulas
    ratio_A_V_D = sympy.simplify(A_D / V_D)
    
    # Define the expected ratio
    expected_ratio = mu_0 * epsilon_0 * v_vec
    
    if ratio_A_V_D != expected_ratio:
        return (f"Incorrect: The formulas in Option D do not satisfy the required physical "
                f"relationship A = (v/c^2)V. The ratio A/V evaluates to {ratio_A_V_D}, "
                f"but it should be mu_0*epsilon_0*v_vec.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_lienard_wiechert_potentials()
print(result)