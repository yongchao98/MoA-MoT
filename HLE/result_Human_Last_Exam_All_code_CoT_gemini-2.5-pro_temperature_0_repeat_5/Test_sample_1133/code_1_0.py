# The problem states that the checkout time for a customer is 10 minutes.
# We are told to assume this is the expected value for the process.
expected_time = 10

# The process described (constant probability of completion per unit time)
# follows an exponential distribution. A key feature of this distribution
# is the "memoryless" property.

# The friend's checkout is a new, independent event from this same distribution.
# Therefore, the expected time for the friend's checkout is equal to the
# mean of the distribution.
expected_time_for_friend = expected_time

# The final equation is simply that the expected time for the friend
# is equal to the given expected time of the process.
print(f"Given Expected Checkout Time E[T] = {expected_time} minutes.")
print(f"The expected time for the friend is also E[T] = {expected_time_for_friend} minutes.")