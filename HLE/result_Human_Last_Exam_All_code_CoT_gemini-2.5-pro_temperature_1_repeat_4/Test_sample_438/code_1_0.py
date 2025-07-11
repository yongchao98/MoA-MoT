def calculate_disease_risk_change():
    """
    This function models the change in disease risk based on antigen presentation probability.
    """
    # A hypothetical baseline probability of presenting a specific self-antigen
    # that leads to an autoimmune response. Let's assume it's very low.
    baseline_presentation_probability = 0.0005

    # The fold increase in presentation probability due to the HLA variant.
    fold_increase = 1000

    # Calculate the new probability of presentation.
    new_presentation_probability = baseline_presentation_probability * fold_increase

    print("### Biological Rationale ###")
    print("Human Leukocyte Antigen (HLA) molecules present protein fragments (antigens) to T-cells.")
    print("An autoimmune disease occurs when the immune system mistakenly attacks the body's own cells after being presented with a 'self-antigen'.")
    print("If an HLA variant makes the presentation of a disease-causing self-antigen more efficient, it amplifies the signal for this mistaken attack.")
    print("\n### Risk Calculation Model ###")
    print(f"Let's assume a baseline probability of this self-antigen being presented is: {baseline_presentation_probability}")
    print(f"The specific HLA variant increases this probability by a factor of: {fold_increase}")
    
    # Per instructions, printing each number in the final equation.
    print("\nThe new, amplified probability is calculated as:")
    print(f"New Probability = Baseline Probability * Fold Increase")
    print(f"{new_presentation_probability:.2f} = {baseline_presentation_probability} * {fold_increase}")
    
    print("\nConclusion: Such a dramatic amplification in presenting a disease-causing self-antigen would almost certainly lead to a stronger, more sustained autoimmune response.")
    print("Therefore, the person's risk of developing the disease would likely increase.")

# Run the calculation and print the results.
calculate_disease_risk_change()