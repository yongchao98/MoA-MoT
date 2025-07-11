def solve_checkout_problem():
    """
    This function solves the supermarket checkout expected time problem.
    """
    # The problem describes a process that follows an exponential distribution.
    # We are given the expected value (mean) of this distribution.
    expected_time_per_customer = 10  # minutes

    # The key insight is the "memoryless" property of the exponential distribution.
    # This means that the time it takes for one customer to check out is
    # completely independent of the time it took for any previous customer.

    # Your friend is starting their checkout process now. Their expected checkout
    # time is a new, independent sample from this distribution.
    # Therefore, their expected time is the same as the overall average.
    expected_time_for_friend = expected_time_per_customer

    # The final equation demonstrates this conclusion.
    # We are asked to find the expected time for the friend (E_friend).
    # We are given the expected time for any customer (E_customer).
    # Because the events are independent: E_friend = E_customer.
    print("The expected checkout time for any customer is given.")
    print(f"E_customer = {expected_time_per_customer} minutes")
    print("\nSince each customer's checkout is an independent event, the expected time for the friend is the same.")
    print("Final Equation:")
    print(f"E_friend = {expected_time_for_friend}")
    print(f"\nThe expected time until your friend is done is {expected_time_for_friend} minutes.")

solve_checkout_problem()