def solve_checkout_problem():
    """
    This function solves the supermarket checkout expected time problem.
    """
    # The problem states that the "10 minutes" figure represents the expected
    # value of the checkout time for a customer.
    # Let E_customer be the expected checkout time for any given customer.
    expected_time_any_customer = 10  # minutes

    # Your friend is about to start their checkout process.
    # Their checkout is a new, independent event that follows the same
    # probability distribution as any other customer's checkout.
    # Therefore, the expected time for your friend is the same as the
    # expected time for any other customer.
    expected_time_friend = expected_time_any_customer

    # The final equation is a simple statement of this fact.
    # Let E_friend be the expected time for the friend.
    print("Let E_friend be the expected checkout time for the friend.")
    print("Let E_customer be the expected checkout time for any customer.")
    print(f"From the problem description, we know: E_customer = {expected_time_any_customer} minutes.")
    print("\nSince the friend's checkout is a new and independent event, it has the same expected value.")
    print("Final Equation: E_friend = E_customer")
    print(f"Result: E_friend = {expected_time_friend} minutes.")

solve_checkout_problem()