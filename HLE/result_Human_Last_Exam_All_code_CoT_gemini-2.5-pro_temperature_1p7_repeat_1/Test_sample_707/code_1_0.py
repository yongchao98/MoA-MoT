def flawed_risk_calculation():
    """
    This function demonstrates Agent A's flawed calculation.
    He incorrectly assumes that the total probability of failure
    is the sum of the probability increases for each quarter.
    """
    num_quarters = 20
    increase_per_quarter = 0.01

    # Create a list of strings representing the increase for each quarter
    equation_terms = [str(increase_per_quarter)] * num_quarters

    # Join the terms with ' + ' to form the equation string
    equation_str = " + ".join(equation_terms)

    # Calculate the flawed result
    flawed_total_probability = increase_per_quarter * num_quarters

    # Print the full, flawed equation and its result
    print("Agent A's flawed calculation was based on the false premise that quarterly risk increases can be linearly added:")
    print(f"{equation_str} = {flawed_total_probability:.1f}")

flawed_risk_calculation()