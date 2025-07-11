import sys

def solve_checkout_problem():
    """
    Calculates the expected checkout time for a friend based on the provided information.
    """
    # The problem states that the average checkout time, based on a prior customer, is 10 minutes.
    # This value represents the expected value of the checkout time distribution.
    expected_time_any_customer = 10  # minutes

    # The friend's checkout is a new, independent event from the same process.
    # The underlying process has a constant probability of completion per second, which
    # implies an exponential distribution. The key property is that each customer's
    # checkout time is drawn from the same distribution.
    # Therefore, the expected checkout time for the friend is the same as for any other customer.
    expected_time_friend = expected_time_any_customer

    print("The expected checkout time for any given customer is 10 minutes.")
    print("Since your friend's checkout is a new instance of the same process, their expected checkout time is also 10 minutes.")
    print("\nThe final equation is:")
    print(f"Expected_Time_Friend = Expected_Time_Any_Customer")
    # Using the print function to display the numbers in the final equation.
    sys.stdout.write(f"{expected_time_friend} = {expected_time_any_customer}\n")

solve_checkout_problem()