def calculate_risk_change():
    """
    This function models how a change in antigen presentation
    can affect a relative risk factor for a disease.
    """
    # Let's assign a baseline relative risk factor for the disease with a common HLA variant.
    # We can set it to 1 for this demonstration.
    baseline_risk_factor = 1

    # The new HLA variant increases the probability of presenting the key self-antigen.
    increase_factor = 1000

    # The new relative risk is the baseline risk multiplied by the increase factor.
    # This is a simplified model to illustrate the concept.
    new_risk_factor = baseline_risk_factor * increase_factor

    print("--- Autoimmune Disease Risk Model ---")
    print("An autoimmune disease can be triggered when the immune system attacks a 'self-antigen'.")
    print("The risk is related to how well an individual's HLA molecules present this self-antigen.")
    print("\nLet's model the change in risk:")
    print(f"Baseline relative risk factor with a common HLA variant: {baseline_risk_factor}")
    print(f"The new variant increases presentation by a factor of: {increase_factor}")
    print("\nFinal Equation for the new relative risk factor:")
    
    # Printing each number in the final equation as requested.
    print(f"{baseline_risk_factor} (baseline) * {increase_factor} (increase factor) = {new_risk_factor} (new risk factor)")
    
    print("\nConclusion: A 1000-fold increase in presenting a disease-causing self-antigen dramatically increases the relative risk of developing the disease.")

calculate_risk_change()