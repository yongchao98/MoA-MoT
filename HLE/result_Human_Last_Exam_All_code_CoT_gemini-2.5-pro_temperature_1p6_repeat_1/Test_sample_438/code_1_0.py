def analyze_risk_change():
    """
    Analyzes and models the change in disease risk due to an HLA variant.
    """
    # 1. Define the variables based on the problem.
    # Let's assume a hypothetical baseline risk for the disease. For example, 0.05%
    baseline_risk = 0.0005

    # The problem states the probability of presentation increases by a factor of 1000.
    fold_increase = 1000

    # 2. Model the new risk.
    # We can create a simple model where the disease risk is directly proportional to
    # the presentation efficiency of the disease-causing antigen.
    new_risk = baseline_risk * fold_increase

    # 3. Display the results and the underlying equation.
    print("Modeling the change in disease risk:")
    print("-" * 35)
    print(f"Assumed Baseline Risk: {baseline_risk} (or {baseline_risk:.2%})")
    print(f"Factor of Increased Presentation: {fold_increase}")

    print("\nTo calculate the new risk, we use the equation:")
    print("New Risk = Baseline Risk * Factor of Increase")
    
    # Per the instructions, printing the equation with the actual numbers.
    print(f"New Risk = {baseline_risk} * {fold_increase}")

    print(f"\nCalculated New Risk: {new_risk} (or {new_risk:.2%})")
    print("-" * 35)

    # 4. State the conclusion.
    if new_risk > baseline_risk:
        print("\nConclusion: The new risk is substantially higher than the baseline risk.")
        print("This means the person's risk of developing the disease would likely INCREASE.")
    else:
        print("\nConclusion: The risk would not increase.")

# Run the analysis
analyze_risk_change()
