import math

def check_answer():
    """
    Checks the correctness of the formula for the maximum speed of a 1D
    relativistic harmonic oscillator. The proposed correct formula is B.

    B) v_max = c * sqrt(1 - 1 / (1 + (k * A^2) / (2 * m * c^2))^2)
    """

    # This is the formula from option B, which is being checked.
    def formula_to_check(m, k, A, c):
        if m <= 0 or c <= 0 or k < 0:
            # Handle invalid physical parameters
            return float('nan')
        try:
            # The term gamma = 1 + (k*A^2)/(2*m*c^2)
            gamma = 1 + (k * A**2) / (2 * m * c**2)
            # The argument of the square root is 1 - 1/gamma^2
            # Since gamma >= 1, the argument is always non-negative.
            v_max_sq_over_c_sq = 1 - 1 / (gamma**2)
            return c * math.sqrt(v_max_sq_over_c_sq)
        except (ValueError, ZeroDivisionError):
            # Catch potential math errors, though the formula is robust
            return float('nan')

    # Test 1: Numerical check against first principles
    m, k, A, c = 1.0, 1.0, 1.0, 2.0
    # From first principles (conservation of energy):
    # E_total = mc^2 + 1/2 * kA^2 = 1*(2^2) + 0.5*1*(1^2) = 4.5
    # E_total = gamma_max * mc^2  =>  gamma_max = E_total / (mc^2) = 4.5 / (1 * 4) = 9/8
    # From gamma_max = 1/sqrt(1 - v_max^2/c^2), we solve for v_max:
    # v_max = c * sqrt(1 - 1/gamma_max^2)
    expected_v_max = c * math.sqrt(1 - 1 / (9/8)**2)
    calculated_v_max = formula_to_check(m, k, A, c)

    if not math.isclose(expected_v_max, calculated_v_max, rel_tol=1e-12):
        return (f"Incorrect. The formula failed the numerical test case derived from energy conservation.\n"
                f"For m={m}, k={k}, A={A}, c={c}:\n"
                f"Expected v_max (from first principles): {expected_v_max}\n"
                f"Calculated v_max (from the given formula): {calculated_v_max}")

    # Test 2: Edge Case - Zero Amplitude
    m, k, A, c = 1.0, 100.0, 0.0, 3e8
    expected_v_max = 0.0
    calculated_v_max = formula_to_check(m, k, A, c)
    if not math.isclose(expected_v_max, calculated_v_max, rel_tol=1e-12):
        return (f"Incorrect. The formula failed the zero amplitude (A=0) edge case.\n"
                f"Expected v_max: {expected_v_max}\n"
                f"Calculated v_max: {calculated_v_max}")

    # Test 3: Edge Case - Zero Spring Constant
    m, k, A, c = 1.0, 0.0, 5.0, 3e8
    expected_v_max = 0.0
    calculated_v_max = formula_to_check(m, k, A, c)
    if not math.isclose(expected_v_max, calculated_v_max, rel_tol=1e-12):
        return (f"Incorrect. The formula failed the zero spring constant (k=0) edge case.\n"
                f"Expected v_max: {expected_v_max}\n"
                f"Calculated v_max: {calculated_v_max}")

    # Test 4: Non-relativistic (classical) limit
    m, k, A, c = 1.0, 100.0, 0.001, 3e8  # Case where v << c
    # Classical formula: 1/2*m*v^2 = 1/2*k*A^2 => v = sqrt(k*A^2/m)
    classical_v_max = math.sqrt(k * A**2 / m)
    relativistic_v_max = formula_to_check(m, k, A, c)
    # For low energies, the relativistic result should be very close to the classical one.
    if not math.isclose(classical_v_max, relativistic_v_max, rel_tol=1e-9):
        return (f"Incorrect. The formula does not correctly reduce to the classical limit for low energies.\n"
                f"Expected classical v_max: {classical_v_max}\n"
                f"Calculated relativistic v_max: {relativistic_v_max}")

    # If all tests pass, the formula is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)