import re

def check_maxwell_monopole_answer():
    """
    Checks the correctness of the LLM's answer about Maxwell's equations and magnetic monopoles.
    """
    # The final answer provided by the LLM
    llm_answer_text = "<<<C>>>"
    
    # --- Problem Constraints and Physics Principles ---
    
    # Principle 1: Gauss's Law for Magnetism (Divergence of B)
    # Standard form (∇ ⋅ B = 0) states no magnetic monopoles exist.
    # If monopoles exist, they are sources of the B field, so ∇ ⋅ B ≠ 0.
    # Conclusion: The law for the divergence of the magnetic field must change.
    divergence_B_changes = True

    # Principle 2: Faraday's Law of Induction (Circulation/Curl of E)
    # By symmetry, if moving electric charges (current) create a circulating B-field (Ampere's Law),
    # then moving magnetic charges (magnetic current) must create a circulating E-field.
    # This requires adding a magnetic current term to Faraday's Law.
    # Conclusion: The law for the circulation of the electric field must change.
    circulation_E_changes = True

    # Principle 3: Gauss's Law for Electricity (Flux/Divergence of E)
    # This law relates the E-field to electric charges. It is unaffected by magnetic charges.
    flux_E_changes = False

    # Principle 4: Ampere-Maxwell Law (Circulation/Curl of B)
    # This law relates the B-field to electric currents and changing E-fields. It is unaffected by magnetic charges.
    circulation_B_changes = False

    # The correct state of changes
    correct_changes = {
        "divergence_of_magnetic_field": divergence_B_changes,
        "circulation_of_electric_field": circulation_E_changes,
        "flux_of_electric_field": flux_E_changes, # Also divergence of E
        "circulation_of_magnetic_field": circulation_B_changes # Also curl of B
    }

    # --- Options as defined in the prompt ---
    options = {
        "A": {
            "description": "The one related to the divergence of the magnetic field.",
            "changes": {
                "divergence_of_magnetic_field": True,
                "circulation_of_electric_field": False,
                "flux_of_electric_field": False,
                "circulation_of_magnetic_field": False
            }
        },
        "B": {
            "description": "The ones related to the divergence and the curl of the magnetic field.",
            "changes": {
                "divergence_of_magnetic_field": True,
                "circulation_of_electric_field": False,
                "flux_of_electric_field": False,
                "circulation_of_magnetic_field": True
            }
        },
        "C": {
            "description": "The ones related to the circulation of the electric field and the divergence of the magnetic field.",
            "changes": {
                "divergence_of_magnetic_field": True,
                "circulation_of_electric_field": True,
                "flux_of_electric_field": False,
                "circulation_of_magnetic_field": False
            }
        },
        "D": {
            "description": "The one related to the circulation of the magnetic field and the flux of the electric field.",
            "changes": {
                "divergence_of_magnetic_field": False,
                "circulation_of_electric_field": False,
                "flux_of_electric_field": True,
                "circulation_of_magnetic_field": True
            }
        }
    }

    # --- Evaluation ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. Expected '<<<X>>>' where X is A, B, C, or D."

    llm_choice = match.group(1)
    
    chosen_option_changes = options[llm_choice]["changes"]

    if chosen_option_changes == correct_changes:
        return "Correct"
    else:
        # Find the discrepancies
        errors = []
        if chosen_option_changes["divergence_of_magnetic_field"] != correct_changes["divergence_of_magnetic_field"]:
            errors.append("The answer incorrectly states whether the law for the 'divergence of the magnetic field' changes. It should change.")
        if chosen_option_changes["circulation_of_electric_field"] != correct_changes["circulation_of_electric_field"]:
            errors.append("The answer incorrectly states whether the law for the 'circulation of the electric field' changes. It should change.")
        if chosen_option_changes["flux_of_electric_field"] != correct_changes["flux_of_electric_field"]:
            errors.append("The answer incorrectly states whether the law for the 'flux of the electric field' changes. It should not change.")
        if chosen_option_changes["circulation_of_magnetic_field"] != correct_changes["circulation_of_magnetic_field"]:
            errors.append("The answer incorrectly states whether the law for the 'circulation of the magnetic field' changes. It should not change.")
        
        # Check for incompleteness
        if not errors: # This happens if the choice is a subset of the correct answer
            if correct_changes["divergence_of_magnetic_field"] and not chosen_option_changes["divergence_of_magnetic_field"]:
                 errors.append("The answer is incomplete. It fails to mention that the law for the 'divergence of the magnetic field' changes.")
            if correct_changes["circulation_of_electric_field"] and not chosen_option_changes["circulation_of_electric_field"]:
                 errors.append("The answer is incomplete. It fails to mention that the law for the 'circulation of the electric field' changes.")

        return f"Incorrect. The chosen option '{llm_choice}' is wrong for the following reason(s):\n- " + "\n- ".join(errors)

# Run the check
result = check_maxwell_monopole_answer()
print(result)