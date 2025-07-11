def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the
    memoryless property of the checkout process.
    """
    # The problem states that the expected checkout time for a customer is 10 minutes.
    # This value is given based on observing a previous customer.
    expected_time_for_any_customer = 10  # minutes

    print("This problem involves a process with a constant probability of completion per unit of time.")
    print("This type of process is known as 'memoryless'.")
    print("\nKey facts:")
    print("1. Each customer's checkout is an independent event.")
    print("2. The process for each customer follows the same probability distribution.")
    print(f"3. We are given that the expected time for this process is {expected_time_for_any_customer} minutes.")

    # The question asks for the expected time for your friend, who is just starting.
    # Since the friend's checkout is a new, independent instance of the same process,
    # its expected value is the same.
    expected_time_for_friend = expected_time_for_any_customer

    print("\nYour friend is a new customer, and their checkout is independent of the previous one.")
    print("Therefore, their expected checkout time is the same as any other customer's.")

    # Print the final equation and answer
    print("\nFinal Calculation:")
    print(f"Expected_Time(Friend) = Expected_Time(Any Customer) = {expected_time_for_friend} minutes")

solve_checkout_problem()