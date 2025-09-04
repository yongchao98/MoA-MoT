import re

def check_answer():
    """
    Checks the correctness of the LLM's answer based on the physics of magnetic monopoles.
    """
    # The final answer provided by the LLM
    llm_answer_key = "C"

    # The options from the question
    options = {
        "A": "The ones related to the divergence and the curl of the magnetic field.",
        "B": "The one related to the divergence of the magnetic field.",
        "C": "The ones related to the circulation of the electric field and the divergence of the magnetic field.",
        "D": "The one related to the circulation of the magnetic field and the flux of the electric field."
    }

    if llm_answer_key not in options:
        return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(options.keys())}."

    answer_text = options[llm_answer_key]

    # --- Define the physical constraints ---

    # Constraint 1: The answer MUST mention the change to the divergence of the magnetic field.
    # This is because ∇ ⋅ B = 0 becomes ∇ ⋅ B = ρ_m.
    must_contain = [
        "divergence of the magnetic field",
        "circulation of the electric field" # This is because ∇ × E = -∂B/∂t becomes ∇ × E = -∂B/∂t - J_m
    ]

    # Constraint 2: The answer MUST NOT mention changes to the other two laws.
    # ∇ ⋅ E = ρ_e / ε₀ and ∇ × B = μ₀(J_e + ε₀ ∂E/∂t) are unaffected.
    must_not_contain = [
        "curl of the magnetic field",       # Part of Ampere-Maxwell Law
        "circulation of the magnetic field", # Same as curl of B
        "flux of the electric field"        # Part of Gauss's Law for Electricity
    ]

    # --- Check the constraints ---

    # Check for required terms
    missing_terms = []
    for term in must_contain:
        # Use regex for word boundary matching to avoid partial matches
        if not re.search(r'\b' + re.escape(term) + r'\b', answer_text):
            missing_terms.append(term)

    if missing_terms:
        return f"Incorrect. The answer is incomplete. It fails to mention that the law related to the '{', '.join(missing_terms)}' would change."

    # Check for forbidden terms
    extra_terms = []
    for term in must_not_contain:
        if re.search(r'\b' + re.escape(term) + r'\b', answer_text):
            extra_terms.append(term)

    if extra_terms:
        return f"Incorrect. The answer incorrectly claims that the law related to the '{', '.join(extra_terms)}' would change, which is not the case."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)