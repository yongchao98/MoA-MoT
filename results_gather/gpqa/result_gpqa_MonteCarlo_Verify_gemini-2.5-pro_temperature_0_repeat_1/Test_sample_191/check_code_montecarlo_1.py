import sympy as sp

def check_correctness_of_e_field_answer(answer_choice: str):
    """
    Checks the correctness of the answer to the electrostatics problem.

    The function works by applying the principles of electrostatic shielding and Gauss's Law
    to derive the correct expression for the electric field outside the conductor. It then
    compares this derived expression with the one corresponding to the provided answer choice.

    Args:
        answer_choice (str): The letter corresponding to the selected answer (e.g., 'A', 'B', 'C', 'D').

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Step 1: Define symbolic variables for the physical quantities ---
    # This allows us to work with the formulas directly.
    q = sp.Symbol('q', positive=True)        # Charge inside the cavity
    epsilon_o = sp.Symbol('epsilon_o', positive=True) # Permittivity of free space
    L = sp.Symbol('L', positive=True)        # Distance from conductor's center to point P
    l = sp.Symbol('l', positive=True)        # Distance from cavity's center to point P
    s = sp.Symbol('s', positive=True)        # Distance between centers
    theta = sp.Symbol('theta')               # Angle

    # Define the electrostatic constant k for brevity
    k = 1 / (4 * sp.pi * epsilon_o)

    # --- Step 2: Apply Physics Principles to derive the correct expression ---
    # Principle 1: Charge Induction.
    # The charge +q inside the cavity induces an equal and opposite charge, -q, on the
    # inner surface (the wall of the cavity).
    # Since the conductor is overall electrically neutral, a charge of +q must appear
    # on the outer surface of the conductor to maintain neutrality.
    charge_on_outer_surface = q

    # Principle 2: Electrostatic Shielding (Faraday Cage).
    # A conducting shell shields its exterior from the electric fields of charges inside it.
    # Therefore, the electric field at point P outside the conductor is produced *only*
    # by the charge distributed on the outer surface of the conductor. The original charge +q
    # and the induced charge -q on the inner surface have no net effect on the field at P.

    # Principle 3: Gauss's Law and Spherical Symmetry (Shell Theorem).
    # The induced charge +q on the outer surface of the *spherical* conductor will
    # distribute itself uniformly over that surface.
    # The electric field outside a uniformly charged spherical shell is identical to the
    # field of a point charge (with the same total charge) located at the center of the sphere.

    # Derivation:
    # - The effective charge creating the external field is `charge_on_outer_surface` = q.
    # - The effective location of this charge is the center of the spherical conductor.
    # - The distance from this center to the point P is L.
    # - Therefore, the magnitude of the electric field E at P is given by Coulomb's Law for a point charge.
    correct_expression = k * charge_on_outer_surface / L**2

    # --- Step 3: Define the expressions for all possible answers ---
    options = {
        'A': k * q / l**2,
        'B': k * q / (l + s * sp.cos(theta))**2,
        'C': k * q / L**2,
        'D': k * q / (l - s * sp.cos(theta))**2
    }

    # --- Step 4: Compare the LLM's answer with the correct expression ---
    if answer_choice not in options:
        return f"Invalid answer choice '{answer_choice}'. Please choose from A, B, C, or D."

    llm_expression = options[answer_choice]

    # The sp.Eq() function performs a structural comparison of the symbolic expressions.
    if sp.Eq(llm_expression, correct_expression):
        return "Correct"
    else:
        reason = (
            "The answer is incorrect.\n"
            "Reasoning: The problem relies on the principle of electrostatic shielding. The conductor isolates the external region from the charge inside the cavity.\n"
            "1. The charge +q in the cavity induces a charge of +q on the *outer surface* of the spherical conductor.\n"
            "2. Due to the spherical symmetry of the conductor, this outer charge +q distributes itself uniformly.\n"
            "3. According to the Shell Theorem (a result of Gauss's Law), the electric field outside a uniformly charged sphere is the same as if all the charge were concentrated in a point charge at its center.\n"
            "4. The center of the conductor is at a distance L from the point P.\n"
            f"Therefore, the correct magnitude of the electric field is E = (1 / (4 * pi * epsilon_o)) * q / L^2. The chosen answer corresponds to option {answer_choice}, which is incorrect."
        )
        return reason

# The provided answer from the other LLM is 'C'.
llm_answer = "C"

# Run the check.
result = check_correctness_of_e_field_answer(llm_answer)

# The final output will be the result of this check.
# If the result is "Correct", we return that.
# If it's an error message, we return the message.
print(result)