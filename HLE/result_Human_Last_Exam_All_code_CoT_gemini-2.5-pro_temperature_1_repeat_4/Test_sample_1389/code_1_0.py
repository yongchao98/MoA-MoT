def analyze_vaccine_efficacy():
    """
    Analyzes if 1-IRR correctly estimates per-exposure efficacy for an all-or-nothing vaccine.
    """
    # Step 1: Define parameters for the "all-or-nothing" vaccine model.
    # Let's assume a vaccine that gives 80% of recipients perfect immunity.
    # The true per-exposure efficacy (VE_p) is the proportion of the vaccinated
    # population that is fully protected.
    v_true_efficacy = 0.8

    # Let's assume an underlying infection hazard rate for any susceptible individual.
    # The exact value doesn't matter as it will cancel out, but we use one for a clear demonstration.
    lambda_hazard = 0.05 # e.g., 5 cases per 100 person-years at risk

    print("--- Step 1: Model Setup ---")
    print(f"The true per-exposure vaccine efficacy (v) is the proportion of vaccinated people who are fully protected.")
    print(f"True VE (v) = {v_true_efficacy}")
    print(f"The underlying hazard of infection (λ) for a susceptible person is {lambda_hazard}.")
    print("-" * 60)

    # Step 2: Calculate the Incidence Rate in the Unvaccinated Group (IR_u).
    # Since everyone in this group is susceptible, their incidence rate is simply lambda.
    IR_u = lambda_hazard
    print("--- Step 2: Calculate Incidence Rate in Unvaccinated Group ---")
    print("The incidence rate in the unvaccinated group (IR_u) equals the base hazard λ.")
    print(f"IR_u = λ = {IR_u}")
    print("-" * 60)

    # Step 3: Calculate the Incidence Rate in the Vaccinated Group (IR_v).
    # This is a weighted average. A proportion 'v' is immune (rate=0) and '1-v' is susceptible (rate=λ).
    IR_v = (v_true_efficacy * 0) + ((1 - v_true_efficacy) * lambda_hazard)
    print("--- Step 3: Calculate Incidence Rate in Vaccinated Group ---")
    print("The rate (IR_v) is a weighted average: (v * rate_protected) + ((1-v) * rate_unprotected)")
    print(f"IR_v = ({v_true_efficacy} * 0) + ((1 - {v_true_efficacy}) * {lambda_hazard})")
    print(f"IR_v = {IR_v}")
    print("-" * 60)

    # Step 4: Calculate the Estimated VE from the Incidence Rate Ratio (IRR).
    IRR = IR_v / IR_u
    VE_estimated = 1 - IRR
    print("--- Step 4: Calculate Estimated VE using 1 - IRR ---")
    print(f"First, calculate the Incidence Rate Ratio (IRR) = IR_v / IR_u")
    print(f"IRR = {IR_v} / {IR_u}")
    print(f"IRR = {IRR}")
    print(f"\nNext, calculate the estimated VE = 1 - IRR")
    print(f"Estimated VE = 1 - {IRR}")
    print(f"Estimated VE = {VE_estimated:.2f}")
    print("-" * 60)

    # Step 5: Compare the true VE with the estimated VE.
    print("--- Step 5: Final Comparison and Conclusion ---")
    print(f"True per-exposure VE (v)      = {v_true_efficacy}")
    print(f"Estimated VE from 1 - IRR     = {VE_estimated:.2f}")

    if abs(VE_estimated - v_true_efficacy) < 1e-9:
        conclusion = "correctly estimates"
    elif VE_estimated > v_true_efficacy:
        conclusion = "overestimates"
    else:
        conclusion = "underestimates"

    print(f"\nConclusion: For an all-or-nothing vaccine, 1 - Incidence Rate Ratio {conclusion} the per-exposure vaccine efficacy.")

analyze_vaccine_efficacy()