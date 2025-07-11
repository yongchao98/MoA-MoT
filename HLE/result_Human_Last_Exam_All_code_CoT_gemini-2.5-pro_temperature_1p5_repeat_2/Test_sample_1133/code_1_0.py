# The problem describes a process with a constant probability of completion per second,
# which corresponds to an exponential distribution. A key property of this distribution
# is that it is "memoryless".

# We are given that the expected checkout time for a customer is 10 minutes.
expected_time_per_customer = 10

# Your friend's checkout is a new, independent event that follows the same
# probability distribution. The time taken by the previous customer sets our
# expectation for the process, but it doesn't influence the outcome for the next person.
# Therefore, the expected checkout time for your friend is the same as the
# established expected time for any customer.

expected_time_friend = expected_time_per_customer

# The final equation is simply stating this equality.
print(f"The final equation is: Expected Time (Friend) = Expected Time (Any Customer)")
print(f"Substituting the given value: Expected Time (Friend) = {expected_time_friend}")
print(f"\nThe expected time until your friend is done is {expected_time_friend} minutes.")
