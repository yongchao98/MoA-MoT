def flawed_risk_calculation():
    """
    This function demonstrates Agent A's flawed reasoning by linearly adding
    quarterly probability increases to calculate the total probability of failure.
    """
    num_quarters = 20
    # A assumed a constant increase of 1% (0.01) in failure probability each quarter.
    quarterly_increase = 0.01

    total_failure_prob = 0
    
    print("Agent A's flawed calculation:")
    
    # We construct the string representing the sum.
    # Note: A's error was to sum the assumed increases, leading to exactly 0.2.
    equation_parts = [str(quarterly_increase) for _ in range(num_quarters)]
    equation_str = " + ".join(equation_parts)
    
    total_failure_prob = num_quarters * quarterly_increase
    
    print(f"P(Failure) = {equation_str} = {total_failure_prob:.2f}")

    print("\nExplanation of the Error:")
    print("Agent A's premise was that the total probability of failure is the sum of the increases in probability from each quarter.")
    print("This is incorrect because probabilities of non-mutually exclusive events cannot be simply added together.")
    print("This method overestimates the risk and is a fundamental misapplication of probability theory.")

flawed_risk_calculation()