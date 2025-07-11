def solve_checkout_problem():
    """
    Calculates the expected checkout time for the friend based on the
    memoryless property of the underlying exponential distribution.
    """
    # The problem states that the expected checkout time for a customer is 10 minutes.
    # This is the mean of the underlying probability distribution.
    # Let's denote the expected time as E.
    E = 10

    # Your friend has just arrived at the cashier, so their process is starting from t=0.
    # For a process described by an exponential distribution (implied by the constant
    # probability of finishing), the expected time from the start is simply the mean of that distribution.
    # This is a direct application of the memoryless property.
    expected_time_for_friend = E

    # Print the explanation and the final equation.
    print(f"The expected checkout time for any customer is given as E = {E} minutes.")
    print("The checkout process is memoryless, meaning the expected future time is independent of the past.")
    print("Since your friend is just starting, their expected checkout time is equal to E.")
    print(f"Final Answer: The expected time for your friend to finish is {expected_time_for_friend} minutes.")

solve_checkout_problem()