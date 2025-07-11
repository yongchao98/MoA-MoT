def solve_checkout_problem():
    """
    Solves the supermarket checkout expected time problem.
    """
    # The problem describes a process where the probability of an event (finishing checkout)
    # is constant over time. This implies the time taken for the event follows an
    # exponential distribution.

    # We are given that the expected checkout time for a customer is 10 minutes.
    # This is the mean of our distribution.
    expected_time_per_customer = 10

    # Each customer's checkout is an independent and identically distributed (i.i.d.) event.
    # This means the time taken by the first customer does not affect the time taken by your friend.
    # Your friend's checkout is a new, independent draw from the same distribution.

    # The question asks for the expected time until your friend is done, from the moment
    # they arrive at the cashier.
    expected_time_for_friend = expected_time_per_customer

    # The final equation is simply stating that the friend's expected time is equal to the
    # known expected time for any customer.
    print(f"The expected checkout time for any given customer is {expected_time_per_customer} minutes.")
    print("Each customer's checkout time is an independent event.")
    print("Therefore, the expected checkout time for your friend is the same.")
    print(f"Final Equation: Expected Time for Friend = {expected_time_for_friend}")
    print(f"\nThe expected time until your friend is done is {expected_time_for_friend} minutes.")

solve_checkout_problem()
<<<C>>>