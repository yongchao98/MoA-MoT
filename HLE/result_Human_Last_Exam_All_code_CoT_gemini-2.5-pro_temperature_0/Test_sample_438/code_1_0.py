def calculate_disease_risk_change():
    """
    This function models the change in probability of presenting a disease-causing
    self-antigen due to an HLA variant.
    """
    # Let's assume the baseline probability of presenting the specific self-antigen
    # in an individual without the variant is very low, for example, 1 in 1,000,000.
    baseline_probability = 0.000001

    # The problem states the variant increases this probability by 1000 fold.
    fold_increase = 1000

    # Calculate the new, higher probability for an individual with the variant.
    variant_probability = baseline_probability * fold_increase

    print("Answering the question: This would likely INCREASE a person's risk of developing the disease.")
    print("A higher probability of presenting a disease-causing self-antigen increases the chance of triggering an autoimmune response.")
    print("\n--- Illustrating the Change in Probability ---")
    print(f"Baseline Probability (without variant): {baseline_probability}")
    print(f"Fold Increase (due to variant): {fold_increase}")
    print("\nFinal Equation:")
    print(f"New Probability = Baseline Probability * Fold Increase")
    # The f-string below prints each number in the final equation.
    print(f"{variant_probability:.4f} = {baseline_probability} * {fold_increase}")
    print(f"\nThe probability increases from 1 in a million to 1 in a thousand, a significant rise in risk.")

calculate_disease_risk_change()