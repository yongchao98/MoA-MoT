def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the
    properties of the exponential distribution.
    """
    # The problem states the checkout process can be modeled by a distribution
    # where the completion probability is constant per unit of time. This implies
    # an exponential distribution for the checkout duration.

    # We are given that the expected value (average time) for a customer's
    # checkout is 10 minutes. This defines the mean of our distribution.
    expected_time_per_customer = 10

    # The friend's checkout is a new, independent event that follows the same
    # probability distribution.
    # Therefore, the expected time for the friend is the same as for any customer.
    expected_time_for_friend = expected_time_per_customer

    print(f"The expected checkout time for any customer is given as: {expected_time_per_customer} minutes.")
    print("Since the friend's checkout is an independent process following the same distribution,")
    print(f"the expected time until the friend is done is also: {expected_time_for_friend} minutes.")

solve_checkout_problem()