def calculate_disease_risk_change():
    """
    This function illustrates how a significant increase in the presentation
    of a disease-causing self-antigen affects disease risk.
    """
    
    # Let's represent the baseline probability of presenting the specific
    # self-antigen as a very small number. This is a hypothetical value for illustration.
    baseline_presentation_probability = 0.000002  # e.g., 2 in a million

    # The HLA variant increases this probability by a given factor.
    fold_increase = 1000

    # Calculate the new probability of presentation with the new HLA variant.
    new_presentation_probability = baseline_presentation_probability * fold_increase
    
    # Explain the biological consequence
    print("The primary question is: Would a 1000-fold increase in presenting a disease-causing self-antigen increase or decrease disease risk?")
    print("\nAnswer: It would significantly INCREASE the risk.\n")
    print("--- Biological Rationale ---")
    print("1. HLA class II molecules present peptides to T-cells to trigger an immune response.")
    print("2. Autoimmune diseases occur when T-cells incorrectly react to the body's own 'self-peptides'.")
    print("3. Increasing the presentation of a 'disease-causing self-peptide' makes it much more likely that self-reactive T-cells will be activated.")
    print("4. This activation leads to an autoimmune attack and, consequently, a higher risk of disease.\n")

    print("--- Illustrative Calculation ---")
    print(f"Let's assume the original probability of presentation is: {baseline_presentation_probability}")
    print(f"The described fold increase is: {fold_increase}")
    
    print("\nThe new, much higher probability is calculated by the equation:")
    # Print each number in the final equation as requested
    print(f"{baseline_presentation_probability} * {fold_increase} = {new_presentation_probability}")
    print("\nThis massive increase in antigen presentation directly translates to a much higher risk of developing the disease.")

# Run the function to display the explanation and calculation.
calculate_disease_risk_change()