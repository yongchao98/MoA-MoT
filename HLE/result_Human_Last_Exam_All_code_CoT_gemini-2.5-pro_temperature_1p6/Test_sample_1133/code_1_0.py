def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the provided information.
    """

    # The problem states the checkout process has a constant completion probability per second.
    # This implies the time duration follows an exponential distribution, which is memoryless.
    
    # We are told that the expected checkout time for a customer is 10 minutes.
    # This value is established from observing a previous customer.
    # Let E_checkout be the expected value of the checkout time distribution.
    E_checkout = 10  # minutes

    # Your friend starts their checkout process. This is a new, independent event.
    # It follows the same probability distribution as any other customer's checkout.
    # Therefore, the expected time for your friend's checkout is the same as the
    # general expected checkout time.
    E_friend_checkout = E_checkout

    print("The expected checkout time for any customer is given.")
    print(f"E_checkout = {E_checkout} minutes")
    print("\nYour friend's checkout is a new process that follows the same distribution.")
    print("Therefore, their expected checkout time is the same.")
    print(f"Expected time for friend = E_checkout = {E_friend_checkout} minutes")

solve_checkout_problem()
<<<C>>>