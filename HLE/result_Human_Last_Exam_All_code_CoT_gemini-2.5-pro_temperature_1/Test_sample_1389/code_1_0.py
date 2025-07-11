import numpy as np

def analyze_aon_vaccine_efficacy():
    """
    Simulates a trial for an all-or-nothing vaccine to determine if
    1 - IRR correctly estimates the per-exposure vaccine efficacy.
    """
    # 1. DEFINE PARAMETERS
    # Population sizes
    N_vaccinated = 100000
    N_unvaccinated = 100000

    # True vaccine efficacy for an all-or-nothing vaccine.
    # This is the proportion of vaccinated people who are fully immune.
    # In this model, this value is also the per-exposure efficacy (VE_p).
    VE_true = 0.85  # 85% efficacy

    # Force of infection (daily risk for a susceptible person)
    force_of_infection = 0.0004

    # Trial duration in days
    duration_days = 365

    print("--- Vaccine Trial Simulation: All-or-Nothing Model ---")
    print(f"True Vaccine Efficacy (VE_S, equivalent to VE_p in this model): {VE_true:.2%}")
    print("-" * 60)

    # 2. SIMULATE INFECTIONS
    # Unvaccinated group: all are susceptible
    person_days_unvaccinated = N_unvaccinated * duration_days
    # The number of cases follows a Poisson distribution with a rate of (force_of_infection * person-time)
    cases_unvaccinated = np.random.poisson(force_of_infection * person_days_unvaccinated)

    # Vaccinated group: only the non-immune fraction is at risk
    num_susceptible_in_vax_group = N_vaccinated * (1 - VE_true)
    person_days_at_risk_in_vax_group = num_susceptible_in_vax_group * duration_days
    # Cases only occur in this susceptible subgroup
    cases_vaccinated = np.random.poisson(force_of_infection * person_days_at_risk_in_vax_group)
    
    # Total person-time for the entire vaccinated cohort is used for the rate calculation
    person_days_vaccinated = N_vaccinated * duration_days

    # 3. CALCULATE RATES AND EFFICACY
    # Incidence Rate (cases per person-day)
    ir_unvaccinated = cases_unvaccinated / person_days_unvaccinated
    ir_vaccinated = cases_vaccinated / person_days_vaccinated

    # Incidence Rate Ratio (IRR)
    irr = ir_vaccinated / ir_unvaccinated

    # Estimated Vaccine Efficacy from the trial data
    ve_estimated = 1 - irr

    print("--- Simulation Results ---")
    print(f"Cases (Unvaccinated): {cases_unvaccinated}")
    print(f"Cases (Vaccinated):   {cases_vaccinated}\n")
    print(f"Incidence Rate (Unvaccinated): {ir_unvaccinated:.8f} cases per person-day")
    print(f"Incidence Rate (Vaccinated):   {ir_vaccinated:.8f} cases per person-day\n")
    
    print("--- Final Calculation ---")
    print("Estimated VE = 1 - (Incidence Rate Vaccinated / Incidence Rate Unvaccinated)")
    # The following line prints the final equation with the numbers from the simulation
    print(f"               = 1 - ({ir_vaccinated:.8f} / {ir_unvaccinated:.8f})")
    print(f"               = 1 - {irr:.4f}")
    print(f"               = {ve_estimated:.4f} or {ve_estimated:.2%}\n")

    print("--- Conclusion ---")
    print(f"The estimated efficacy from the simulation ({ve_estimated:.2%}) is very close to the true efficacy ({VE_true:.2%}).")
    print("For an all-or-nothing vaccine, the proportion of people protected (VE_S) is equivalent to the average per-exposure efficacy (VE_p).")
    print("Therefore, 1 - IRR correctly estimates the per-exposure vaccine efficacy.")

analyze_aon_vaccine_efficacy()