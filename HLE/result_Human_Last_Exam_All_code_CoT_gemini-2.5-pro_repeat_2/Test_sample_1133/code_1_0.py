def solve_checkout_problem():
    """
    Calculates the expected checkout time for the friend based on the
    properties of the exponential distribution.
    """
    # The problem states that the expected checkout time for a customer is 10 minutes.
    # This is the mean of the underlying exponential distribution.
    average_checkout_time = 10

    # The friend's checkout is a new, independent event drawn from the same distribution.
    # The exponential distribution is "memoryless", meaning the duration of one customer's
    # checkout does not influence the next customer's checkout time.
    # Therefore, the expected time for the friend is the same as the average time.
    expected_time_for_friend = average_checkout_time

    print("This problem describes a process governed by the exponential distribution.")
    print("A key property of this distribution is that it is memoryless.")
    print("This means each customer's checkout is an independent event.")
    print(f"We are given that the average (expected) checkout time is {average_checkout_time} minutes.")
    print("Therefore, your friend's expected checkout time is also the same.")
    print("\nFinal Equation:")
    print(f"Expected Time for Friend = Average Checkout Time = {expected_time_for_friend}")


solve_checkout_problem()
<<<C>>>