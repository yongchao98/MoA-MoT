import sys

def solve_checkout_problem():
    """
    Solves the supermarket checkout expected time problem.
    """
    # The problem states that the observed checkout time of 10 minutes
    # should be treated as the expected value for the process.
    # Let E(T) be the expected time for any customer to complete checkout.
    expected_time_for_any_customer = 10  # in minutes

    # This type of process, where the probability of an event occurring per unit of time
    # is constant, is described by a distribution with the "memoryless property"
    # (like the Geometric or Exponential distribution).

    # The memoryless property means that the expected future waiting time is
    # independent of the time that has already passed.
    # P(T > t + s | T > s) = P(T > t)

    # Your friend is just beginning their checkout at the cashier.
    # Therefore, we want to find the expected time from the start of their process.
    # Because the process is memoryless, the experience of the previous customer
    # does not influence the expected time for your friend. The expected time for
    # your friend is the same as the fundamental expected time for any customer.

    expected_time_for_friend = expected_time_for_any_customer

    # Print the explanation and the final equation.
    print("This problem is an application of the memoryless property of certain probability distributions.")
    print(f"The given expected checkout time for any customer is {expected_time_for_any_customer} minutes.")
    print("Since your friend is just starting, their expected waiting time is the same as the base case.")
    print("\nFinal Calculation:")
    print(f"Expected time for friend = Expected time for any customer")
    print(f"Expected time for friend = {expected_time_for_friend} minutes")

# Execute the function to print the solution.
# This avoids running the code if the script is imported into another module.
if __name__ == "__main__":
    solve_checkout_problem()
    # The answer corresponds to 10 minutes.
    # Appending the final choice as per the format instruction.
    # We use sys.stdout.write to prevent adding a newline, ensuring the format is exact.
    sys.stdout.write("\n<<<C>>>")