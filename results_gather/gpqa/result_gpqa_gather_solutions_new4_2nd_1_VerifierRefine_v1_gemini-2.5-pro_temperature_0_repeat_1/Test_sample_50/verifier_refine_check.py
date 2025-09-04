import math

def check_correctness():
    """
    This function checks the correctness of the provided answer for the potential energy
    of a point charge near a grounded conducting sphere.

    The function does the following:
    1. Defines Python functions for each of the multiple-choice options (A, B, C, D).
    2. Defines a function for the ground truth formula, which is derived from first principles
       (the method of images).
    3. The final answer from the LLM is identified as 'A'.
    4. It uses a set of sample numerical values for the physical quantities (k, q, R, d)
       that satisfy the problem's constraints (d > R).
    5. It calculates the potential energy using the formula from the chosen answer ('A') and
       the ground truth formula.
    6. It compares the two results. If they match, the answer is correct. Otherwise, it's incorrect.
    """

    # The final answer provided by the LLM analysis.
    llm_final_answer_choice = "A"

    # --- Formula Definitions ---

    # The ground truth is derived from the method of images: U = (1/2) * q * V_induced
    # where V_induced = -k*q*R / (d^2 - R^2)
    # So, the correct formula is: U = - (1/2) * k * q^2 * R / (d^2 - R^2)
    def ground_truth_formula(k, q, R, d):
        # This is the correct formula based on physics principles.
        return -0.5 * k * q**2 * R / (d**2 - R**2)

    # Option A: U = - (1/2) * k * q^2 * R / (d^2 - R^2)
    def formula_A(k, q, R, d):
        return -0.5 * k * q**2 * R / (d**2 - R**2)

    # Option B: U = - (1/2) * k * q^2 * d / (d^2 + R^2)
    def formula_B(k, q, R, d):
        return -0.5 * k * q**2 * d / (d**2 + R**2)

    # Option C: U = - k * q^2 * d / (d^2 - R^2)
    def formula_C(k, q, R, d):
        return -1.0 * k * q**2 * d / (d**2 - R**2)

    # Option D: U = - (1/2) * k * q^2 * R^2 / (d^2 - R^2)
    def formula_D(k, q, R, d):
        return -0.5 * k * q**2 * R**2 / (d**2 - R**2)

    # --- Verification Logic ---

    # Map the answer choices to their respective functions
    formulas = {
        "A": formula_A,
        "B": formula_B,
        "C": formula_C,
        "D": formula_D
    }
    
    chosen_formula_func = formulas.get(llm_final_answer_choice)

    # Choose some sample values for the variables, respecting the constraint d > R.
    k = 1.0  # Using k=1 for simplicity, as it's a common factor.
    q = 2.0
    R = 3.0
    d = 5.0

    # Check the physical constraint
    if not d > R:
        return f"Invalid sample values: The distance 'd' ({d}) must be greater than the radius 'R' ({R})."
    
    # Calculate the result from the chosen formula and the ground truth
    chosen_result = chosen_formula_func(k, q, R, d)
    ground_truth_result = ground_truth_formula(k, q, R, d)

    # Compare the results using math.isclose for robust floating-point comparison
    if math.isclose(chosen_result, ground_truth_result):
        return "Correct"
    else:
        # If incorrect, provide a detailed reason.
        formula_strings = {
            "A": "U = - (1/2) * k * q^2 * R / (d^2 - R^2)",
            "B": "U = - (1/2) * k * q^2 * d / (d^2 + R^2)",
            "C": "U = - k * q^2 * d / (d^2 - R^2)",
            "D": "U = - (1/2) * k * q^2 * R^2 / (d^2 - R^2)"
        }
        
        reason = f"The answer '{llm_final_answer_choice}' is incorrect.\n"
        reason += f"The formula for option {llm_final_answer_choice} is: {formula_strings[llm_final_answer_choice]}\n"
        reason += f"The correct formula, derived from the method of images, is: {formula_strings['A']}\n"
        reason += f"Using sample values (k={k}, q={q}, R={R}, d={d}):\n"
        reason += f"  - The result from the chosen answer's formula ({llm_final_answer_choice}) is: {chosen_result}\n"
        reason += f"  - The result from the correct ground truth formula is: {ground_truth_result}\n"
        reason += "The results do not match."
        return reason

# Execute the check and print the result.
print(check_correctness())