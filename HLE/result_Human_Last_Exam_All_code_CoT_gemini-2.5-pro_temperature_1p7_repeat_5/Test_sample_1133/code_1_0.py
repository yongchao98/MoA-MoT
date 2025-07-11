def solve_checkout_problem():
    """
    Calculates the expected checkout time for your friend based on the problem description.
    """
    
    # The problem states that the expected time for a customer's entire checkout process is 10 minutes.
    # This value sets the mean (expected value) for the underlying probability distribution,
    # which is an exponential distribution.
    expected_time_given = 10  # in minutes

    # Your friend starts their checkout process from the beginning. Their checkout is a new,
    # independent event that follows the same distribution as any other customer.
    # Therefore, their expected time is simply the mean of that distribution.
    expected_time_friend = expected_time_given

    print("The checkout process is memoryless, best described by an exponential distribution.")
    print(f"The given expected value for any customer's full checkout is: {expected_time_given} minutes.")
    print("Your friend is a new customer, so their process starts from the beginning.")
    print("Their expected checkout time is therefore the same as the given expected value.")
    
    # Final equation showing that the friend's expected time equals the given expected time.
    print(f"Final Equation: Expected_Time_Friend ({expected_time_friend}) = Expected_Time_Given ({expected_time_given})")
    
    print(f"\nThe expected time until your friend is done is {expected_time_friend} minutes.")

solve_checkout_problem()
<<<C>>>