import sympy
from sympy.physics.wigner import clebsch_gordan

def check_angular_momentum_probability():
    """
    This function verifies the answer to a quantum mechanics problem involving
    the addition of angular momentum.

    The problem asks for the probability of measuring m1=-1 and m2=-1
    for a system in the state |l1=1, l2=1, l=2, m=-1>.

    The probability is given by the square of the Clebsch-Gordan coefficient:
    P = |<l1, m1; l2, m2 | l, m>|^2
    """
    try:
        # --- Define quantum numbers from the problem ---

        # Initial coupled state: |l1, l2, l, m>
        l1 = 1
        l2 = 1
        l = 2
        m = -1

        # Desired measurement outcome (final uncoupled state): |l1, m1; l2, m2>
        # The eigenvalues of L_z are m*hbar.
        # An eigenvalue of -hbar corresponds to a quantum number of -1.
        m1 = -1
        m2 = -1

        # The provided answer is A, which corresponds to a probability of 0.
        expected_probability = 0

        # --- Calculate the correct probability ---

        # A fundamental selection rule for Clebsch-Gordan coefficients is that
        # they are non-zero only if m = m1 + m2. This reflects the conservation
        # of the z-component of angular momentum.
        if m != m1 + m2:
            # If the rule is violated, the coefficient is zero, and so is the probability.
            calculated_probability = 0
        else:
            # If the rule is satisfied, we would calculate the coefficient explicitly.
            # This part is not strictly necessary for this problem but makes the checker robust.
            cg_coeff = clebsch_gordan(l1, l2, l, m1, m2, m)
            # The probability is the square of the coefficient.
            calculated_probability = float(cg_coeff**2)

        # --- Compare and return the result ---

        # Use a small tolerance for floating-point comparison.
        if abs(calculated_probability - expected_probability) < 1e-9:
            return "Correct"
        else:
            reason = (
                f"The provided answer corresponds to a probability of {expected_probability}, "
                f"but the calculated probability is {calculated_probability}.\n"
                f"The reason is the violation of a fundamental selection rule.\n"
                f"The total z-component of angular momentum for the initial state is m = {m}.\n"
                f"The desired measurement outcome has individual components m1 = {m1} and m2 = {m2}, "
                f"which sum to m1 + m2 = {m1 + m2}.\n"
                f"For the probability to be non-zero, it must be that m = m1 + m2. "
                f"Since {m} != {m1 + m2}, the probability of this outcome is exactly 0."
            )
            return reason

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check and print the result.
print(check_angular_momentum_probability())