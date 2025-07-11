import numpy as np

def analyze_vaccine_efficacy():
    """
    Analyzes if 1 - Incidence Rate Ratio (IRR) correctly estimates the
    efficacy of an "all-or-nothing" vaccine.
    """

    # --- 1. Define model parameters ---
    # This is the base infection rate in the general population.
    # Let's assume 5 out of 100 unvaccinated people get the disease per year.
    lambda_rate = 0.05

    # This is the true efficacy of the vaccine.
    # For an all-or-nothing vaccine, this means 80% of vaccinated people
    # are fully protected, and 20% are not protected at all.
    ve_true = 0.80

    print("This simulation checks if 1 - IRR is a good estimator for an 'all-or-nothing' vaccine.")
    print("\n--- Model Setup ---")
    print(f"True Vaccine Efficacy (VE_true): {ve_true:.2f} (This means {ve_true*100:.0f}% of vaccine recipients are fully immune)")
    print(f"Base infection rate in unvaccinated population (lambda): {lambda_rate:.2f}")

    # --- 2. Calculate theoretical incidence rates ---

    # In the unvaccinated group, the incidence rate is the base rate.
    ir_unvaccinated = lambda_rate
    print("\nStep 1: Calculate the Incidence Rate in the Unvaccinated Group (IR_u).")
    print(f"IR_u = lambda = {ir_unvaccinated:.4f}")

    # In the vaccinated group, only the non-protected portion (1 - ve_true) is at risk.
    ir_vaccinated = lambda_rate * (1 - ve_true)
    print("\nStep 2: Calculate the Incidence Rate in the Vaccinated Group (IR_v).")
    print("This rate is the base rate times the fraction of people not protected by the vaccine.")
    print(f"IR_v = lambda * (1 - VE_true) = {lambda_rate:.2f} * (1 - {ve_true:.2f}) = {ir_vaccinated:.4f}")

    # --- 3. Calculate the Incidence Rate Ratio (IRR) ---
    irr = ir_vaccinated / ir_unvaccinated
    print("\nStep 3: Calculate the Incidence Rate Ratio (IRR).")
    print(f"IRR = IR_v / IR_u = {ir_vaccinated:.4f} / {ir_unvaccinated:.4f} = {irr:.2f}")

    # --- 4. Calculate the vaccine efficacy from the IRR ---
    ve_measured = 1 - irr
    print("\nStep 4: Calculate the measured vaccine efficacy (VE_measured) from the IRR.")
    print(f"VE_measured = 1 - IRR = 1 - {irr:.2f} = {ve_measured:.2f}")

    # --- 5. Compare true vs. measured efficacy ---
    print("\n--- Conclusion ---")
    print(f"The true vaccine efficacy we started with was: {ve_true:.2f}")
    print(f"The efficacy calculated using 1 - IRR is: {ve_measured:.2f}")

    if np.isclose(ve_true, ve_measured):
        print("\nFinal Answer: 1 - IRR correctly estimates the per-exposure vaccine efficacy.")
    elif ve_measured > ve_true:
        print("\nFinal Answer: 1 - IRR overestimates the per-exposure vaccine efficacy.")
    else:
        print("\nFinal Answer: 1 - IRR underestimates the per-exposure vaccine efficacy.")


if __name__ == "__main__":
    analyze_vaccine_efficacy()
