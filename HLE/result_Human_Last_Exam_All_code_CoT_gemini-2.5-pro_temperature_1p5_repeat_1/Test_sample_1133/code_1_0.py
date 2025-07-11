def solve_supermarket_puzzle():
    """
    Solves the supermarket checkout expected time puzzle.
    """
    # Step 1: Identify the given expected value for a single checkout process.
    # The problem states that the expected checkout time for a customer is 10 minutes.
    # This information establishes the parameter of the underlying probability distribution.
    # The phrase "constant probability ... per second" implies an exponential distribution.
    # For an exponential distribution, the mean (expected value) is the key parameter.
    expected_time_any_customer = 10  # in minutes

    # Step 2: Analyze the situation for the friend.
    # The friend is starting a new, independent checkout process.
    # This new process follows the same exponential distribution as any other customer.
    # Therefore, the expected time for the friend's checkout is the mean of that distribution.
    expected_time_friend = expected_time_any_customer

    # Step 3: Print the reasoning and the final equation.
    print("The key to this problem lies in understanding the properties of the checkout process.")
    print("1. The process has a 'constant probability of completion per second', which means the time duration follows an exponential distribution.")
    print("2. We are given that the expected value (mean) of this distribution is 10 minutes.")
    print("3. Your friend is starting a new, independent checkout. Their checkout time is a new random sample from this same distribution.")
    print("4. Therefore, the expected time for your friend is simply the mean of the distribution.")
    print("\nFinal Calculation:")
    print(f"Expected Time for Friend = Expected Time for Any Customer = {expected_time_friend} minutes")

solve_supermarket_puzzle()