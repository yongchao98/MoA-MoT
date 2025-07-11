def solve_checkout_time():
    """
    Calculates the expected checkout time for the friend based on the problem's description.
    """
    # The problem states that the expected checkout time for a customer is 10 minutes.
    # This value is given directly.
    # In an exponential distribution, E[T] = 1/lambda. We are given E[T].
    expected_time_customer = 10

    # The friend's checkout process is an independent event that follows the same
    # probability distribution as any other customer.
    # Therefore, the expected time for the friend's checkout is the same.
    expected_time_friend = expected_time_customer

    print("The problem describes a process where the duration follows an exponential distribution.")
    print(f"The given expected checkout time for any single customer is: {expected_time_customer} minutes.")
    print("The friend's checkout is a new, independent instance of this process.")
    print("\nTherefore, the expected time until the friend is done is given by the equation:")
    print(f"Expected Time (Friend) = Expected Time (Any Customer)")
    print(f"Expected Time (Friend) = {expected_time_friend}")

solve_checkout_time()