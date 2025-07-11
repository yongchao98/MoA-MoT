def analyze_disease_risk():
    """
    Analyzes how a change in antigen presentation probability affects disease risk.
    """
    # Step 1: Define the initial parameters based on the problem.
    # We'll use a baseline probability of 1 (as an arbitrary unit) for simplicity.
    baseline_presentation_probability = 1
    
    # The variant increases this probability by a factor of 1000.
    fold_increase = 1000
    
    # Step 2: Calculate the new probability.
    new_presentation_probability = baseline_presentation_probability * fold_increase
    
    # Step 3: Explain the logic and show the calculation.
    print("Analyzing the effect of the HLA variant on disease risk:")
    print("-" * 50)
    print("The risk of certain autoimmune diseases is linked to how strongly the immune system responds to a 'self-antigen'.")
    print("If an HLA variant presents this self-antigen more effectively, it increases the chance of an autoimmune reaction.")
    print("\nLet's quantify the change:")
    
    # As requested, outputting each number in the final equation.
    print(f"Baseline Presentation Probability (arbitrary units): {baseline_presentation_probability}")
    print(f"Fold Increase from HLA Variant: {fold_increase}")
    print(f"Resulting New Presentation Probability: {baseline_presentation_probability} * {fold_increase} = {new_presentation_probability}")
    print("-" * 50)
    
    # Step 4: Draw a conclusion based on the result.
    if new_presentation_probability > baseline_presentation_probability:
        conclusion = "INCREASE"
    else:
        conclusion = "DECREASE or stay the same"
        
    print(f"Conclusion: Since the probability of presenting the disease-causing self-antigen is {fold_increase} times higher, it would likely {conclusion} a person's risk of developing the disease.")

# Run the analysis
analyze_disease_risk()