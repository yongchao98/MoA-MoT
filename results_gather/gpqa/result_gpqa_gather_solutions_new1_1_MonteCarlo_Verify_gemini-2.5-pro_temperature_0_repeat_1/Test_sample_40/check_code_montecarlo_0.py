import math

def check_answer_correctness():
    """
    Checks the correctness of the proposed answer for the speed of light in a moving medium.

    The question asks for the speed of light in a glass (index n) moving at velocity v.
    The speed of light in vacuum c=1.
    The proposed answer is A: (1 + n*v) / (n + v)
    """

    # The formula from the proposed answer 'A'
    formula_A = lambda n, v: (1 + n * v) / (n + v)

    # --- Test 1: Limiting Case for Vacuum (n=1) ---
    # If the medium is a vacuum, the speed of light must be c=1 for any observer.
    n_vacuum = 1.0
    v_test = 0.5  # An arbitrary velocity < c
    expected_speed_in_vacuum = 1.0
    calculated_speed_in_vacuum = formula_A(n_vacuum, v_test)

    if not math.isclose(calculated_speed_in_vacuum, expected_speed_in_vacuum):
        return (f"Incorrect. The answer fails the vacuum constraint (n=1).\n"
                f"For n={n_vacuum}, the speed should be {expected_speed_in_vacuum}, "
                f"but the formula gives {calculated_speed_in_vacuum:.4f}.")

    # --- Test 2: Limiting Case for Stationary Medium (v=0) ---
    # If the medium is stationary, the speed of light should be c/n = 1/n.
    n_glass = 1.5  # A typical refractive index for glass
    v_stationary = 0.0
    expected_speed_in_stationary_medium = 1.0 / n_glass
    calculated_speed_in_stationary_medium = formula_A(n_glass, v_stationary)

    if not math.isclose(calculated_speed_in_stationary_medium, expected_speed_in_stationary_medium):
        return (f"Incorrect. The answer fails the stationary medium constraint (v=0).\n"
                f"For n={n_glass}, the speed should be 1/n = {expected_speed_in_stationary_medium:.4f}, "
                f"but the formula gives {calculated_speed_in_stationary_medium:.4f}.")

    # --- Final Check ---
    # The formula from answer A is identical to the one derived from the relativistic
    # velocity addition formula: u = (1 + n*v) / (n + v).
    # Since it passes both critical physical constraints, it is correct.
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)