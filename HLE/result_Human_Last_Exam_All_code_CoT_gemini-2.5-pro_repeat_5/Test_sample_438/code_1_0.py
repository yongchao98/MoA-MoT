def calculate_risk_change(baseline_presentation, increase_factor):
    """
    Calculates and explains the effect of increased antigen presentation on disease risk.

    Args:
        baseline_presentation (float): The initial relative probability of presenting a self-antigen.
        increase_factor (int): The fold increase in presentation probability due to the HLA variant.
    """

    # Calculate the new presentation probability
    new_presentation = baseline_presentation * increase_factor

    print("Understanding the Impact of Increased Antigen Presentation:")
    print("-" * 55)
    print("1. HLA molecules present peptides (antigens) to T-cells.")
    print("2. The immune system is normally tolerant to 'self-antigens'.")
    print("3. An autoimmune disease can occur if this tolerance breaks and the immune system attacks 'self'.")
    print("\nScenario:")
    print(f"An HLA variant increases the presentation of a disease-causing self-antigen by {increase_factor}-fold.")

    print("\nIllustrative Calculation:")
    print("Let's model the relative change in presentation probability.")
    # The prompt requires printing each number in the final equation.
    print(f"New Presentation Probability = Baseline Probability * Increase Factor")
    print(f"New Presentation Probability = {baseline_presentation} * {increase_factor} = {new_presentation}")

    print("\nConclusion:")
    print(f"The new relative probability of presentation is {new_presentation}, which is {increase_factor} times higher than the baseline.")
    print("This dramatically increases the chance of activating self-reactive T-cells, which would likely INCREASE a person's risk for the disease.")

# --- Parameters ---
# Let's set a baseline relative probability to 1 for simplicity.
baseline_prob = 1
# The fold increase given in the problem.
factor = 1000

# Run the calculation and explanation
calculate_risk_change(baseline_prob, factor)