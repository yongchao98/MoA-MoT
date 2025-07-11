def calculate_risk_change(baseline_risk, fold_increase):
    """
    Calculates and displays the effect of a fold increase on a baseline value.

    Args:
        baseline_risk (float): The initial baseline risk, represented as a unit.
        fold_increase (float): The factor by which the risk increases.
    """
    # Calculate the new risk after the increase
    new_risk = baseline_risk * fold_increase

    # The prompt requires outputting each number in the final equation.
    # We will print the equation to show the relationship between the numbers.
    print(f"If we model the baseline probability of presentation as: {baseline_risk}")
    print(f"And the fold increase caused by the genetic variant is: {fold_increase}")
    print("The new relative probability is the product of these two numbers.")
    print(f"Final Equation: {baseline_risk} * {fold_increase} = {new_risk}")
    print(f"\nThis {fold_increase}-fold increase in presenting a disease-causing self-antigen would likely increase disease risk.")

# Define the baseline and the fold increase from the problem description
baseline_probability_unit = 1
increase_factor = 1000

# Run the function to display the calculation
calculate_risk_change(baseline_probability_unit, increase_factor)
