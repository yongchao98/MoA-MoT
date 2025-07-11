import sys

def calculate_risk_change():
    """
    This function calculates and explains the impact of increased antigen presentation on disease risk.
    """
    # Let's assign a baseline arbitrary unit for the probability of presentation.
    baseline_presentation_probability = 1

    # The variant increases this probability by a given factor.
    fold_increase = 1000

    # Calculate the new probability relative to the baseline.
    new_presentation_probability = baseline_presentation_probability * fold_increase

    print(f"If the original probability of presenting a specific self-antigen is: {baseline_presentation_probability}")
    print(f"And a genetic variant causes a fold increase of: {fold_increase}")
    # We must explicitly print the numbers in the final equation as requested.
    print(f"The new relative probability is the original value multiplied by the fold increase:")
    sys.stdout.write(f"{baseline_presentation_probability} * {fold_increase} = {new_presentation_probability}\n")


    print("\nConclusion:")
    print("This dramatically higher level of self-antigen presentation makes it much more likely")
    print("that an autoimmune response will be triggered against the body's own tissues.")
    print("Therefore, this would likely INCREASE a person's risk of developing the disease.")

calculate_risk_change()