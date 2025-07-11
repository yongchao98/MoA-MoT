def explain_rmbs_risk_factors():
    """
    This function conceptually models how originator quality acts as a root cause
    for the decline in non-agency RMBS value during the 2004-2008 period.
    """

    # --- Input Factors ---
    # Let's represent originator quality on a scale of 1 (Poor) to 10 (Excellent).
    # During the 2004-2008 crisis, many originators had very poor quality.
    originator_quality_score = 2  # Represents a poor quality originator

    # --- How Originator Quality Influences Key Risk Metrics ---
    # A poor originator (low score) leads to higher expected defaults and lower recoveries.
    # We model this with simple formulas for illustration.
    # A high-quality originator might see a 2% default rate, a poor one might see 30%+.
    expected_default_rate = 0.5 / originator_quality_score  # Results in 25%

    # Poor originators approved high loan-to-value loans, leading to low recoveries in foreclosure.
    # A high-quality originator might recover 70%, a poor one might recover 40% or less.
    expected_recovery_rate = 0.2 * originator_quality_score # Results in 40%

    # --- Calculating the RMBS Value Conceptually ---
    face_value = 100.0
    loss_given_default = 1.0 - expected_recovery_rate
    total_expected_loss = expected_default_rate * loss_given_default
    final_rmbs_value = face_value * (1.0 - total_expected_loss)

    # --- Explanation ---
    print("### Conceptual Model: Originator Quality as the Root Cause ###\n")
    print(f"The core problem in the 2004-2008 crisis was the quality of the underlying loans, which was determined by the loan originator.")
    print(f"Let's assume a conceptual 'Originator Quality Score' of {originator_quality_score} (out of 10).\n")

    print("Step 1: Originator quality determines the Expected Default Rate.")
    print(f"A low score leads to lax underwriting and thus a high default rate.")
    print(f"Formula: 0.5 / Originator Quality Score = Expected Default Rate")
    print(f"Calculation: 0.5 / {originator_quality_score} = {expected_default_rate:.2f} or {expected_default_rate:.0%}\n")

    print("Step 2: Originator quality determines the Expected Recovery Rate.")
    print(f"Lax underwriting (e.g., high loan-to-value) means less money is recovered after a default.")
    print(f"Formula: 0.2 * Originator Quality Score = Expected Recovery Rate")
    print(f"Calculation: 0.2 * {originator_quality_score} = {expected_recovery_rate:.2f} or {expected_recovery_rate:.0%}\n")

    print("Step 3: Calculate the final estimated value of the RMBS.")
    print("The value is the Face Value minus the Total Expected Loss.")
    print("Equation: Value = Face_Value * (1 - (Expected_Default_Rate * (1 - Expected_Recovery_Rate)))")
    print(f"Final Calculation: {face_value:.2f} * (1 - ({expected_default_rate:.2f} * (1 - {expected_recovery_rate:.2f}))) = {final_rmbs_value:.2f}\n")

    print("Conclusion: As shown, the 'Originator Quality' is the root factor that drives the critical variables (defaults, recoveries) which ultimately determine the RMBS value.")
    print("Therefore, the quality of the loan issuer and RMBS originator was the most likely root cause factor.")


explain_rmbs_risk_factors()
<<<F>>>