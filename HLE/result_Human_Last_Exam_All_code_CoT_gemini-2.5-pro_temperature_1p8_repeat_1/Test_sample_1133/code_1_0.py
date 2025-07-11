def solve_checkout_problem():
    """
    Solves the checkout line problem using the properties of the exponential distribution.

    The problem states that the probability of completing checkout is constant per second.
    This describes a process modeled by the exponential distribution, which has a key
    property: it is "memoryless".

    The memoryless property means that the expected future waiting time is independent
    of the time already spent.

    We are given that the expected checkout time for any customer is 10 minutes. This
    sets the mean of our distribution.
    """

    # The given expected value for any customer's checkout process.
    # This is the mean (μ) of the underlying exponential distribution.
    given_expected_time = 10  # in minutes

    # Your friend's checkout is a new, independent event from the same distribution.
    # Therefore, the expected time for your friend to check out is the same as the
    # established mean of the distribution.
    expected_time_for_friend = given_expected_time

    # Print the explanation and the final equation.
    print("The checkout process is memoryless, modeled by an exponential distribution.")
    print("The given expected time for any customer defines the mean of this distribution.")
    print(f"Given Expected Time = {given_expected_time} minutes")
    print("\nSince the friend's checkout is a new, independent event from the same distribution, their expected time is the same.")
    print("\nFinal Equation:")
    print(f"Expected Time for Friend = Given Expected Time")
    print(f"Expected Time for Friend = {expected_time_for_friend} minutes")

solve_checkout_problem()