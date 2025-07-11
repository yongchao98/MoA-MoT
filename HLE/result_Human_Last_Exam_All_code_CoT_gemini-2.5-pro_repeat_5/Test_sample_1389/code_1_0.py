import sys
# Redirect print to a string buffer to avoid asking user to copy paste
import io
old_stdout = sys.stdout
sys.stdout = mystdout = io.StringIO()

def analyze_all_or_nothing_vaccine():
    """
    Analyzes if 1-IRR correctly estimates per-exposure VE for an all-or-nothing vaccine.
    """
    # --- 1. Define Parameters ---
    # True proportion of vaccinated individuals who are fully protected.
    # This is the true "per-exposure efficacy" (VE_p) at the population level.
    ve_aon = 0.8  # 80% efficacy

    # Hazard rate (force of infection) for any susceptible individual (cases per person-year).
    lambda_susceptible = 0.05

    # Population sizes and time period
    n_vaccinated = 100000
    n_unvaccinated = 100000
    time_years = 1.0

    print("--- Model Parameters ---")
    print(f"True All-or-Nothing Efficacy (VE_p): {ve_aon:.0%}")
    print(f"Hazard Rate for Susceptibles: {lambda_susceptible}")
    print(f"Population Size (Vaccinated & Unvaccinated): {n_vaccinated}")
    print("-" * 26)

    # --- 2. Calculate Expected Cases and Incidence Rates ---

    # Unvaccinated Group: All are susceptible.
    person_years_unvaccinated = n_unvaccinated * time_years
    cases_unvaccinated = n_unvaccinated * lambda_susceptible * time_years
    ir_unvaccinated = cases_unvaccinated / person_years_unvaccinated

    # Vaccinated Group:
    # A fraction (ve_aon) is fully protected (0 cases).
    # The remaining fraction (1 - ve_aon) is fully susceptible.
    n_vaccinated_susceptible = n_vaccinated * (1 - ve_aon)
    person_years_vaccinated = n_vaccinated * time_years
    cases_vaccinated = n_vaccinated_susceptible * lambda_susceptible * time_years
    ir_vaccinated = cases_vaccinated / person_years_vaccinated

    print("\n--- Calculations ---")
    print(f"Incidence Rate in Unvaccinated (IR_u): {ir_unvaccinated:.4f} cases per person-year")
    print(f"Incidence Rate in Vaccinated (IR_v): {ir_vaccinated:.4f} cases per person-year")

    # --- 3. Calculate Vaccine Efficacy from IRR ---
    if ir_unvaccinated > 0:
        irr = ir_vaccinated / ir_unvaccinated
        ve_from_irr = 1 - irr
    else:
        ve_from_irr = None # Cannot calculate if there's no disease in unvaccinated group

    print(f"\nIncidence Rate Ratio (IRR = IR_v / IR_u): {irr:.2f}")
    print(f"Calculated VE (1 - IRR): {ve_from_irr:.2f}")

    # --- 4. Compare and Conclude ---
    print("\n--- Conclusion ---")
    print(f"The true per-exposure efficacy (VE_p) is {ve_aon:.2f}.")
    print(f"The efficacy calculated as 1 - IRR is {ve_from_irr:.2f}.")

    # The final equation as requested.
    print("\nFinal Equation:")
    print(f"1 - (IR_v / IR_u) = 1 - ({ir_vaccinated:.4f} / {ir_unvaccinated:.4f}) = {ve_from_irr:.2f}")

    if round(ve_from_irr, 4) == round(ve_aon, 4):
        print("\nResult: 1 - IRR correctly estimates the per-exposure vaccine efficacy.")
    elif ve_from_irr > ve_aon:
        print("\nResult: 1 - IRR overestimates the per-exposure vaccine efficacy.")
    else:
        print("\nResult: 1 - IRR underestimates the per-exposure vaccine efficacy.")

analyze_all_or_nothing_vaccine()

# Restore standard output and print the captured string
sys.stdout = old_stdout
print(mystdout.getvalue())