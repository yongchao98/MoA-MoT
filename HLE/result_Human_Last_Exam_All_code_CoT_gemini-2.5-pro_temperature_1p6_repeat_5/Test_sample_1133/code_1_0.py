def solve_checkout_problem():
    """
    Solves the checkout line probability puzzle.
    """

    # Step 1: Define the given information from the problem.
    # We are told that the expected checkout time for a customer is 10 minutes.
    # This value represents the mean of the underlying probability distribution.
    expected_time_for_any_customer = 10  # in minutes

    # Step 2: Explain the underlying principle.
    # The problem describes a process where the probability of completion is constant per unit of time.
    # This is a classic signature of the Exponential Distribution.
    # For an exponential distribution, the expected value (mean) is constant for any new, independent event.
    print("The checkout process is modeled by an exponential distribution because the completion probability is constant per second.")

    # Step 3: Apply the principle to the friend's situation.
    # The friend is starting a new checkout process. This is an independent event.
    # Therefore, the expected time for the friend's checkout is simply the mean of the distribution.
    expected_time_for_friend = expected_time_for_any_customer
    print(f"The expected checkout time for any given customer is {expected_time_for_any_customer} minutes.")
    print("Your friend's checkout is a new, independent event that follows the same distribution.")

    # Step 4: Output the final answer and the simple equation.
    print("\nFinal Equation:")
    print(f"Expected Time for Friend = Expected Time for Any Customer = {expected_time_for_friend} minutes")

solve_checkout_problem()
<<<C>>>