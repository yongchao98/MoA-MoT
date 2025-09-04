import sympy

def check_physics_answer():
    """
    Checks the correctness of the answer to the relativistic velocity addition problem.
    """
    # Define symbolic variables for refractive index (n) and velocity (v).
    # c (speed of light in vacuum) is taken as 1.
    n, v = sympy.symbols('n v', real=True, positive=True)
    c = 1

    # --- Step 1: Derive the correct formula from first principles ---
    # The speed of light in the rest frame of the glass (S') is u_prime.
    # By definition, n = c / u_prime, so u_prime = c / n.
    u_prime = c / n

    # The relativistic velocity addition formula for co-linear motion is:
    # u = (v_frame + u_prime) / (1 + (v_frame * u_prime) / c**2)
    # Here, v_frame is the velocity of the glass, which is 'v'.
    # The light beam and glass move in the same direction.
    correct_formula_derived = (v + u_prime) / (1 + (v * u_prime) / c**2)

    # Substitute c=1 and u_prime = 1/n into the derived formula
    correct_formula_derived = correct_formula_derived.subs({c: 1})
    
    # Simplify the derived formula algebraically
    simplified_correct_formula = sympy.simplify(correct_formula_derived)

    # --- Step 2: Define the formulas for the given options ---
    options = {
        'A': (1 + n*v) / (n + v),
        'B': (1 + n*v) / (n - v),
        'C': (1 - n*v) / (n + v),
        'D': sympy.Integer(1)
    }

    # The final answer given in the prompt is 'A'
    chosen_answer_label = 'A'
    chosen_answer_formula = options[chosen_answer_label]

    # --- Step 3: Check if the chosen answer's formula matches the correctly derived formula ---
    # We check if the difference between the two formulas simplifies to zero.
    if sympy.simplify(simplified_correct_formula - chosen_answer_formula) != 0:
        return (f"Incorrect. The formula for the chosen answer {chosen_answer_label}, which is "
                f"{chosen_answer_formula}, does not match the correctly derived formula from "
                f"special relativity, which simplifies to {simplified_correct_formula}.")

    # --- Step 4: Perform sanity checks on the formula ---
    # Sanity Check 1: If the glass is stationary (v=0), the speed should be c/n = 1/n.
    check1_result = sympy.simplify(simplified_correct_formula.subs(v, 0))
    if check1_result != 1/n:
        return (f"Incorrect. The formula {simplified_correct_formula} fails the sanity check for v=0. "
                f"Expected 1/n, but got {check1_result}.")

    # Sanity Check 2: If the medium is a vacuum (n=1), the speed should be c = 1.
    check2_result = sympy.simplify(simplified_correct_formula.subs(n, 1))
    if check2_result != 1:
        return (f"Incorrect. The formula {simplified_correct_formula} fails the sanity check for n=1. "
                f"Expected 1, but got {check2_result}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_physics_answer()
print(result)