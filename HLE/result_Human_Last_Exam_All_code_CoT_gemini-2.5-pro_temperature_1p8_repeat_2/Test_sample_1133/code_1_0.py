def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the problem description.
    """
    # Step 1: The problem states that the expected checkout time for a customer is given.
    # This value defines the mean of the underlying probability distribution for all customers.
    expected_time_for_any_customer = 10  # in minutes

    print("The checkout time for any customer follows a specific probability distribution.")
    print(f"We are given that the expected value (mean) of this distribution is {expected_time_for_any_customer} minutes.")
    print("\nLet E[T_customer] represent this expected time.")
    print(f"E[T_customer] = {expected_time_for_any_customer}")

    # Step 2: The friend is another customer. Their checkout time is an independent
    # event from the same distribution.
    # Therefore, the expected time for the friend is the same as for any other customer.
    expected_time_for_friend = expected_time_for_any_customer

    print("\nYour friend's checkout time is a new, independent sample from the same distribution.")
    print("Therefore, the expected time for your friend is the same as the established mean.")
    print("\nLet E[T_friend] represent the friend's expected checkout time.")
    print("The final equation is: E[T_friend] = E[T_customer]")
    
    # Step 3: Print the final answer.
    print(f"\nFinal Answer: {expected_time_for_friend} minutes")

solve_checkout_problem()
<<<C>>>