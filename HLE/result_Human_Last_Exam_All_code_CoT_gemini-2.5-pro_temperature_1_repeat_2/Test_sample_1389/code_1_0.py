def analyze_all_or_nothing_vaccine():
    """
    Analyzes if 1 - IRR correctly estimates the efficacy of an all-or-nothing vaccine.
    """

    # --- Step 1: Define a hypothetical scenario ---
    # True per-exposure vaccine efficacy (VE_p).
    # This represents the fraction of vaccinated people who are fully protected.
    ve_p = 0.90  # 90% efficacy

    # Rate of exposure to the pathogen (e.g., exposures per year).
    exposure_rate = 5.0

    # Per-exposure infection probability for a susceptible person.
    infection_prob_per_exposure = 0.1

    print("--- Scenario ---")
    print(f"True Per-Exposure Vaccine Efficacy (VE_p): {ve_p:.2f}")
    print(f"Exposure Rate (c): {exposure_rate}")
    print(f"Infection Probability per Exposure (λ): {infection_prob_per_exposure}")
    print("-" * 20)

    # --- Step 2: Calculate Incidence Rates ---
    # Incidence rate in the unvaccinated group (IR_u = c * λ)
    ir_unvaccinated = exposure_rate * infection_prob_per_exposure

    # Incidence rate in the vaccinated group (IR_v = (1 - VE_p) * c * λ)
    # This is because only the fraction (1 - VE_p) is susceptible.
    ir_vaccinated = (1 - ve_p) * exposure_rate * infection_prob_per_exposure

    print("--- Calculations ---")
    print(f"Incidence Rate in Unvaccinated (IR_u): {exposure_rate} * {infection_prob_per_exposure} = {ir_unvaccinated:.4f}")
    print(f"Incidence Rate in Vaccinated (IR_v): (1 - {ve_p}) * {exposure_rate} * {infection_prob_per_exposure} = {ir_vaccinated:.4f}")

    # --- Step 3: Calculate IRR and Estimated Efficacy ---
    # Incidence Rate Ratio (IRR = IR_v / IR_u)
    irr = ir_vaccinated / ir_unvaccinated

    # Estimated Vaccine Efficacy (VE_est = 1 - IRR)
    ve_estimated = 1 - irr

    print(f"\nIncidence Rate Ratio (IRR): {ir_vaccinated:.4f} / {ir_unvaccinated:.4f} = {irr:.4f}")
    print(f"Estimated Efficacy (1 - IRR): 1 - {irr:.4f} = {ve_estimated:.4f}")
    print("-" * 20)

    # --- Step 4: Compare and Conclude ---
    print("--- Conclusion ---")
    print(f"True Efficacy (VE_p): {ve_p:.2f}")
    print(f"Estimated Efficacy (1 - IRR): {ve_estimated:.2f}")

    if abs(ve_p - ve_estimated) < 1e-9:
        print("\nResult: 1 - IRR correctly estimates the per-exposure vaccine efficacy for an all-or-nothing vaccine.")
    elif ve_estimated > ve_p:
        print("\nResult: 1 - IRR overestimates the per-exposure vaccine efficacy for an all-or-nothing vaccine.")
    else:
        print("\nResult: 1 - IRR underestimates the per-exposure vaccine efficacy for an all-or-nothing vaccine.")

if __name__ == "__main__":
    analyze_all_or_nothing_vaccine()