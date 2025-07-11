def analyze_vaccine_efficacy():
    """
    Analyzes if 1 - Incidence Rate Ratio (IRR) correctly estimates
    the efficacy of an 'all-or-nothing' vaccine.
    """
    # Step 1: Define a hypothetical scenario
    # VE_s: True vaccine efficacy (proportion of vaccinated who are fully protected)
    # This is the "per-exposure" efficacy averaged over the entire vaccinated population.
    VE_s = 0.90  # 90% efficacy

    # lambda: Force of infection (rate of new infections in a fully susceptible population)
    lambda_rate = 0.05  # 5 cases per 100 person-years

    print("--- Scenario ---")
    print(f"True Vaccine Efficacy (proportion protected, VE_s): {VE_s:.2f}")
    print(f"Force of Infection (lambda): {lambda_rate:.2f}\n")

    # Step 2: Calculate Incidence Rates
    # Incidence rate in the unvaccinated (everyone is susceptible)
    IR_u = lambda_rate
    
    # Incidence rate in the vaccinated
    # Only the proportion (1 - VE_s) is susceptible
    IR_v = lambda_rate * (1 - VE_s)
    
    print("--- Calculations ---")
    print(f"Incidence Rate in Unvaccinated (IR_u): {IR_u:.4f}")
    print(f"Incidence Rate in Vaccinated (IR_v): {IR_v:.4f}\n")
    
    # Step 3: Calculate Estimated VE
    # Incidence Rate Ratio
    IRR = IR_v / IR_u
    
    # Estimated VE from IRR
    VE_estimated = 1 - IRR
    
    print("--- Final Equation ---")
    # We are calculating: 1 - (IR_v / IR_u)
    print(f"Estimated VE = 1 - ({IR_v:.4f} / {IR_u:.4f})")
    print(f"Estimated VE = 1 - {IRR:.2f}")
    print(f"Estimated VE = {VE_estimated:.2f}\n")
    
    # Step 4: Compare and Conclude
    print("--- Conclusion ---")
    print(f"The true efficacy was {VE_s:.2f}, and the estimated efficacy is {VE_estimated:.2f}.")
    
    if VE_s == VE_estimated:
        print("Therefore, 1 - IRR correctly estimates the vaccine efficacy for an all-or-nothing vaccine.")
    elif VE_s > VE_estimated:
        print("Therefore, 1 - IRR underestimates the vaccine efficacy.")
    else:
        print("Therefore, 1 - IRR overestimates the vaccine efficacy.")

# Run the analysis
analyze_vaccine_efficacy()