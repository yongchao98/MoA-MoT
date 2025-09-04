import math

def check_chemistry_calculation():
    """
    This function checks the correctness of the provided answer for the chemistry problem.
    It recalculates the concentration of orthophosphate ions based on the given data and compares it to the selected option.
    """
    # --- Define constants and initial values from the question ---
    volume_cm3 = 200.00
    mass_kh2po4 = 1.00  # in grams
    mw_kh2po4 = 136.09  # in g/mol
    mass_na2hpo4_2h2o = 1.00  # in grams
    mw_na2hpo4_2h2o = 177.99  # in g/mol
    ka2 = 6.2e-8
    ka3 = 1.8e-12

    # --- Define the options as provided in the question ---
    options = {
        'A': 2.81e-7,
        'B': 5.48e-7,
        'C': 3.97e-7,
        'D': 6.24e-7
    }

    # The final answer provided by the LLM to be checked
    llm_final_answer_key = 'D'

    # --- Perform the calculation step-by-step ---

    # 1. Calculate moles of the buffer components
    moles_h2po4_minus = mass_kh2po4 / mw_kh2po4
    moles_hpo4_2minus = mass_na2hpo4_2h2o / mw_na2hpo4_2h2o

    # 2. Calculate initial concentrations in Liters
    volume_L = volume_cm3 / 1000.0
    conc_h2po4_minus = moles_h2po4_minus / volume_L
    conc_hpo4_2minus = moles_hpo4_2minus / volume_L

    # 3. Determine the hydrogen ion concentration [H+] using the Ka2 equilibrium
    # Equilibrium: H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻
    # Ka₂ = [H⁺][HPO₄²⁻] / [H₂PO₄⁻]
    # [H⁺] = Ka₂ * ([H₂PO₄⁻] / [HPO₄²⁻])
    conc_h_plus = ka2 * (conc_h2po4_minus / conc_hpo4_2minus)

    # 4. Calculate the orthophosphate ion [PO₄³⁻] concentration using the Ka3 equilibrium
    # Equilibrium: HPO₄²⁻ ⇌ H⁺ + PO₄³⁻
    # Ka₃ = [H⁺][PO₄³⁻] / [HPO₄²⁻]
    # [PO₄³⁻] = Ka₃ * ([HPO₄²⁻] / [H⁺])
    calculated_conc_po4_3minus = ka3 * (conc_hpo4_2minus / conc_h_plus)

    # --- Verify the correctness of the LLM's answer ---
    
    # Get the numerical value corresponding to the LLM's chosen answer key
    llm_answer_value = options.get(llm_final_answer_key)

    if llm_answer_value is None:
        return f"The provided answer key '{llm_final_answer_key}' is not a valid option (A, B, C, or D)."

    # Check if the calculated result is close to the value of the chosen option.
    # A relative tolerance of 1% is used to account for potential rounding in the options.
    if math.isclose(calculated_conc_po4_3minus, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # If not correct, find the closest option to the calculated value for a more detailed reason.
        closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_conc_po4_3minus))
        
        reason = (
            f"The provided answer '{llm_final_answer_key}' is incorrect. "
            f"The step-by-step calculation yields an orthophosphate concentration of {calculated_conc_po4_3minus:.3e} M. "
            f"This calculated value is closest to option {closest_option_key} ({options[closest_option_key]:.3e} M), "
            f"not option {llm_final_answer_key} ({llm_answer_value:.3e} M)."
        )
        return reason

# Run the check and print the result
print(check_chemistry_calculation())