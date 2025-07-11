def assess_disease_risk():
    """
    Models the change in disease risk based on an increased probability of
    presenting a disease-causing self-antigen.
    """
    # A hypothetical baseline risk for the disease. Let's assume it's 0.05%
    baseline_risk = 0.0005

    # The factor by which the HLA variant increases the probability of presentation.
    fold_increase = 1000

    # Calculate the new risk. We assume the risk is directly proportional
    # to the probability of antigen presentation.
    new_risk = baseline_risk * fold_increase

    print("Analyzing the impact of an HLA variant on disease risk...")
    print(f"Hypothetical Baseline Risk: {baseline_risk}")
    print(f"Fold Increase in Antigen Presentation: {fold_increase}")
    print("-" * 30)
    
    # Per the instructions, we output the numbers used in the final equation.
    print(f"The new risk is calculated by multiplying the baseline risk by the fold increase.")
    print(f"Equation: New Risk = Baseline Risk * Fold Increase")
    print(f"Calculation: {new_risk:.2f} = {baseline_risk} * {fold_increase}")
    print("-" * 30)


    # Compare the new risk to the baseline risk to draw a conclusion.
    if new_risk > baseline_risk:
        print("Conclusion: The new risk is substantially higher than the baseline risk.")
        print("This would likely INCREASE a person's risk of developing the disease.")
    elif new_risk < baseline_risk:
        print("Conclusion: The new risk is lower than the baseline risk.")
        print("This would likely DECREASE a person's risk of developing the disease.")
    else:
        print("Conclusion: The risk remains unchanged.")

assess_disease_risk()