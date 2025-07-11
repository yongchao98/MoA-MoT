# The problem describes a scenario where the waiting time for a checkout follows an exponential distribution.
# We are given the expected value (mean) of this distribution.

# Expected checkout time for any given customer, as established by the problem.
# This is the mean of the underlying exponential distribution.
mean_checkout_time = 10  # in minutes

# The checkout process for your friend is a new, independent event.
# It is governed by the same probability distribution as any other customer.
# Due to the properties of independent events, the expected checkout time for your friend
# is the same as the mean of the distribution.
expected_time_for_friend = mean_checkout_time

# The final equation shows that the expected time for the friend is equal to the given mean checkout time.
print(f"Expected Time for Friend = {expected_time_for_friend} minutes")