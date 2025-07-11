# The problem states that the checkout time for a customer is 10 minutes.
# It clarifies that this value should be treated as the expected value (mean) of the process.
given_expected_time = 10  # in minutes

# The problem describes a process with a constant probability of completion per unit of time.
# This implies the checkout time follows an exponential distribution, which is "memoryless".
# This means the expected time for any new, independent event (like your friend's checkout)
# is simply the mean of the distribution, regardless of previous events.

# Therefore, the expected time for your friend is equal to the given expected time.
expected_time_for_friend = given_expected_time

# The final equation is a direct assignment based on the properties of the process.
print(f"The expected checkout time for any customer is given as {given_expected_time} minutes.")
print("Since the process is memoryless, the expected time for your friend is the same.")
print(f"Final Equation: Expected Time for Friend = {expected_time_for_friend}")
print(f"The expected time until your friend is done is {expected_time_for_friend} minutes.")
