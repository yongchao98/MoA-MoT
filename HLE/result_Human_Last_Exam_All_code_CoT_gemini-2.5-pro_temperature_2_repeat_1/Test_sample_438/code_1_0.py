def calculate_risk_change(base_risk_percentage, fold_increase):
    """
    Calculates and prints the new disease risk based on a fold increase in
    antigen presentation.

    Args:
        base_risk_percentage (float): The initial risk of the disease as a percentage.
        fold_increase (int): The factor by which antigen presentation increases.
    """
    # Convert percentage to a probability decimal
    initial_risk = base_risk_percentage / 100.0

    # The new risk is directly proportional to the increase in presentation
    new_risk = initial_risk * fold_increase

    print("---Modeling the Change in Disease Risk---")
    print(f"A person's initial risk with a standard HLA gene is {base_risk_percentage}%.")
    print(f"An HLA variant increases the presentation of a disease-causing self-antigen by a factor of {fold_increase}.")
    print("\nAssuming risk is proportional to presentation, we can model the new risk:")
    # Printing the equation with each number as requested
    print(f"New Risk = Initial Risk * Fold Increase")
    print(f"New Risk = {initial_risk} * {fold_increase} = {new_risk}")

    # Convert the new risk back to a percentage for the final statement
    new_risk_percentage = new_risk * 100
    print(f"\nThe new risk is {new_risk_percentage}%, a dramatic increase.")
    print("\nConclusion: This would likely increase a person's risk of developing the disease.")

# ---Parameters---
# Let's assume a hypothetical base risk for the general population.
# For example, Type 1 Diabetes has a lifetime risk of about 0.5% in the general population.
# We'll use this as our example base risk.
base_risk_percentage = 0.5

# The fold increase given in the problem
fold_increase_factor = 1000

# Run the calculation
calculate_risk_change(base_risk_percentage, fold_increase_factor)
