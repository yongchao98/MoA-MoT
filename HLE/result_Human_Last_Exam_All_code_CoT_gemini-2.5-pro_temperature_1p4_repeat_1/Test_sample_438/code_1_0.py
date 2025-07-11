def calculate_disease_risk():
    """
    Calculates and explains the change in disease risk based on a change
    in HLA self-antigen presentation.
    """
    # 1. Define the baseline risk.
    # Let's assume a hypothetical baseline risk of developing the disease
    # for an individual with a common HLA variant.
    baseline_risk = 0.0005

    # 2. Define the risk multiplier from the problem statement.
    # The new variant increases the probability of presentation 1000 fold.
    risk_multiplier = 1000

    # 3. Calculate the new risk.
    new_risk = baseline_risk * risk_multiplier

    # 4. Print the explanation and the calculation.
    print("Thinking Process:")
    print("The presentation of a disease-causing self-antigen by an HLA molecule is a key step in developing an autoimmune disease.")
    print("If a genetic variant makes this presentation much more efficient, the probability of an immune attack against the self-antigen increases.")
    print("\nLet's model this with a simple calculation:")
    print(f"Assume a baseline risk with a normal HLA variant is: {baseline_risk}")
    print(f"The new variant increases this by a factor of: {risk_multiplier}")
    print("\nThe new risk is calculated by multiplying the baseline risk by this factor.")
    
    # As requested, printing the final equation with each number.
    print(f"\nFinal Equation: New Risk = {baseline_risk} * {risk_multiplier}")
    
    print(f"Calculated New Risk = {new_risk}")

    print("\nConclusion:")
    if new_risk > baseline_risk:
        print(f"Since the new risk ({new_risk}) is substantially higher than the baseline risk ({baseline_risk}), this would likely INCREASE the person's risk of developing the disease.")
    else:
        # This case is not expected based on the problem
        print("The risk would not increase.")

# Execute the function
calculate_disease_risk()
