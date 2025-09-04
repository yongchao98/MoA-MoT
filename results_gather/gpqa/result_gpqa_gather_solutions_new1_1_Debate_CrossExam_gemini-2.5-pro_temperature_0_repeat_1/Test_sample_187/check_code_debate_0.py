import math

def check_rhombohedral_distance():
    """
    Calculates the interplanar distance for a rhombohedral crystal and checks
    it against the provided answer.
    """
    # --- Define problem constraints and parameters ---
    # Interatomic distance (lattice parameter) in Angstrom
    a = 10.0
    # Angle between axes in degrees
    alpha_deg = 30.0
    # Miller indices for the plane
    h, k, l = 1, 1, 1

    # The multiple-choice options from the question
    options = {
        'A': 10.05,
        'B': 9.54,
        'C': 8.95,
        'D': 9.08
    }

    # The final answer to be checked
    final_answer_key = 'B'
    final_answer_value = options[final_answer_key]

    # --- Perform the calculation ---
    # For a rhombohedral crystal, the general formula for interplanar distance d_hkl is:
    # 1/d^2 = [ (h^2+k^2+l^2)sin^2(a) + 2(hk+kl+lh)(cos^2(a)-cos(a)) ] / [ a^2(1-3cos^2(a)+2cos^3(a)) ]
    # For the specific (111) plane, this simplifies to a much more stable formula:
    # d_111 = a * sqrt( (1 + 2*cos(a)) / 3 )
    # Using the simplified formula is preferred to avoid floating point errors from the complex one.

    # Convert angle to radians for use in Python's math library
    alpha_rad = math.radians(alpha_deg)

    # Calculate the interplanar distance using the simplified formula
    try:
        calculated_d = a * math.sqrt((1 + 2 * math.cos(alpha_rad)) / 3)
    except ValueError:
        return "Calculation error: The value inside the square root was negative."

    # --- Verify the correctness of the answer ---
    # We check if the calculated value is close to the provided answer's value.
    # A tolerance is used for floating-point comparisons. 0.01 is reasonable
    # given the precision of the options.
    tolerance = 0.01

    if abs(calculated_d - final_answer_value) <= tolerance:
        return "Correct"
    else:
        # If the answer is incorrect, find the correct option and provide a reason.
        correct_key = ''
        min_diff = float('inf')
        for key, value in options.items():
            diff = abs(calculated_d - value)
            if diff < min_diff:
                min_diff = diff
                correct_key = key

        reason = (
            f"Incorrect. The calculated interplanar distance is {calculated_d:.4f} Angstrom. "
            f"The provided answer was {final_answer_value} Angstrom (Option {final_answer_key}). "
            f"The calculated value is closest to Option {correct_key} ({options[correct_key]} Angstrom), "
            f"with a difference of {abs(calculated_d - options[correct_key]):.4f}."
        )
        return reason

# Run the check and print the result
result = check_rhombohedral_distance()
print(result)