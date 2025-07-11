import numpy as np

def analyze_vaccine_efficacy():
    """
    Analyzes and explains if 1-IRR overestimates, underestimates, or correctly
    estimates per-exposure vaccine efficacy for an all-or-nothing vaccine.
    """
    # Step 1: Define parameters for our hypothetical scenario.
    # VE_s is the fraction of vaccinated people who are fully protected (susceptibility is zero).
    # This is the defining parameter of an "all-or-nothing" vaccine.
    ve_s_true = 0.90

    # h is the baseline hazard or incidence rate of infection for any susceptible person.
    # We can use a representative value, e.g., 5 cases per 100 person-years.
    # Let's use cases per person-day for simplicity. h = 0.05 / 365 = 0.000137
    h_unvaccinated = 0.005 # Base incidence rate for susceptible individuals

    # Step 2: Define the 'true' per-exposure vaccine efficacy (VE_p) for this model.
    # For the VE_s fraction, protection is 100%. For the 1-VE_s fraction, it's 0%.
    # The average per-exposure efficacy across the whole group is the weighted average.
    ve_p_true = ve_s_true * 1.0 + (1.0 - ve_s_true) * 0.0

    print("--- Theoretical Framework ---")
    print(f"1. We define an 'all-or-nothing' vaccine with a true efficacy (VE_s) of {ve_s_true:.2f}.")
    print(f"   This means {ve_s_true*100}% of vaccinees are fully immune and the rest have no protection.")
    print(f"2. The 'true' average per-exposure efficacy (VE_p) is the average protection across all vaccinees:")
    print(f"   VE_p = ({ve_s_true:.2f} * 1.0) + ({(1-ve_s_true):.2f} * 0.0) = {ve_p_true:.2f}")

    # Step 3: Calculate the incidence rates for the two groups based on the model.
    # The incidence rate in the unvaccinated (placebo) group is the baseline hazard.
    ir_unvaccinated = h_unvaccinated
    
    # The incidence rate in the vaccinated group is the weighted average hazard.
    # (VE_s * 0% hazard) + ((1-VE_s) * 100% hazard)
    ir_vaccinated = (ve_s_true * 0.0) + ((1.0 - ve_s_true) * h_unvaccinated)
    
    print("\n--- Simulating Study Results ---")
    print(f"3. The incidence rate in the unvaccinated group (IR_u) would be the baseline hazard: {ir_unvaccinated:.4f}")
    print(f"4. The incidence rate in the vaccinated group (IR_v) is the average of the immune and non-immune subgroups:")
    print(f"   IR_v = ({ve_s_true:.2f} * 0.0) + ({(1.0 - ve_s_true):.2f} * {h_unvaccinated:.4f}) = {ir_vaccinated:.4f}")

    # Step 4: Calculate the observed vaccine efficacy from the Incidence Rate Ratio (IRR).
    irr = ir_vaccinated / ir_unvaccinated
    ve_from_irr = 1.0 - irr

    print("\n--- Final Calculation and Conclusion ---")
    print(f"5. From the rates, we calculate the Incidence Rate Ratio (IRR):")
    print(f"   IRR = IR_v / IR_u = {ir_vaccinated:.4f} / {ir_unvaccinated:.4f} = {irr:.2f}")

    print("\n6. The vaccine efficacy (VE) measured in the study is 1 - IRR.")
    # The final equation with numbers is printed here, as requested.
    print(f"   VE = 1 - {irr:.2f} = {ve_from_irr:.2f}")

    print("\nComparison:")
    print(f"Measured VE from study (1 - IRR): {ve_from_irr:.2f}")
    print(f"True per-exposure VE (VE_p):      {ve_p_true:.2f}")

    if np.isclose(ve_from_irr, ve_p_true):
        print("\nResult: The measured VE (1 - IRR) correctly estimates the per-exposure VE.")
    elif ve_from_irr > ve_p_true:
        print("\nResult: The measured VE (1 - IRR) overestimates the per-exposure VE.")
    else:
        print("\nResult: The measured VE (1 - IRR) underestimates the per-exposure VE.")

if __name__ == '__main__':
    analyze_vaccine_efficacy()
<<<C>>>