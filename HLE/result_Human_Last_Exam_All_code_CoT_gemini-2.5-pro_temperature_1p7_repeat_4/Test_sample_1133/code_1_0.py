import sys

def solve_checkout_time():
    """
    This function solves the supermarket checkout probability problem.
    """
    # The problem provides the expected value for a customer's checkout time.
    # This is the mean of the underlying probability distribution for checkout times.
    given_expected_time = 10  # in minutes

    # The friend's checkout is a new, independent event. The problem states
    # the process is the same for everyone. This means the friend's checkout time
    # is drawn from the exact same probability distribution.
    #
    # Therefore, the expected time for the friend's checkout is identical to the
    # expected time for any other customer.

    friend_expected_time = given_expected_time

    # The problem can be represented by a simple identity equation, where the expected
    # time for the friend is equal to the known expected time of the process.
    # We will print out the components of this logical step.

    print("The problem states that the expected checkout time for a customer is 10 minutes.")
    print("Each customer's checkout is an independent event drawn from the same distribution.")
    print("\nLet E_friend be the expected time for the friend.")
    print("Let E_given be the given expected time for any customer.")
    print("\nThe relationship is:")
    print("E_friend = E_given")

    # Displaying the equation with the numbers as requested.
    print("\nSubstituting the values:")
    print(f"{friend_expected_time} = {given_expected_time}")

    # Final conclusion.
    print(f"\nTherefore, the expected time until your friend is done and ready to go is {friend_expected_time} minutes.")

# Execute the solution.
if __name__ == "__main__":
    # The problem is conceptual, but the function below prints the logical steps and answer.
    solve_checkout_time()