def flawed_risk_assessment():
    """
    This function demonstrates Agent A's flawed risk calculation,
    which is based on an incorrect premise about how to cumulate probability.
    """

    # Parameters from Agent A's perspective
    probability_increase_per_quarter = 0.01
    num_years = 5
    quarters_per_year = 4

    # Total number of periods (quarters) over the goal's duration
    total_quarters = num_years * quarters_per_year

    # Agent A's flawed calculation: He incorrectly adds the probability increase
    # for each quarter to arrive at the final cumulative probability.
    final_failure_prob = probability_increase_per_quarter * total_quarters

    print("Agent A's Flawed Calculation")
    print("----------------------------")
    print("The agent's false premise is that cumulative probability increases additively.")
    print("This leads to a simple but incorrect multiplication or summation.")
    print("\nThe flawed equation used is the sum of the probability increases for each quarter:")

    # Build and print the full equation to show every number
    equation_parts = [str(probability_increase_per_quarter)] * total_quarters
    full_equation_str = " + ".join(equation_parts)

    print(f"{full_equation_str} = {final_failure_prob:.2f}")
    print(f"\nThis results in an incorrectly calculated final failure probability of {final_failure_prob:.0%}.")

flawed_risk_assessment()