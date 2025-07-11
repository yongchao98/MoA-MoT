# The problem describes a process with a constant probability of completion per unit of time.
# This implies the checkout time T follows an exponential distribution.

# We are told that the expected value of a customer's checkout time is 10 minutes.
# For an exponential distribution, the expected value (mean) is the parameter that defines it.
expected_checkout_time_for_any_customer = 10

# The problem asks for the expected time for your friend, starting from when they reach the cashier.
# Since the process is the same for everyone, your friend's checkout time
# follows the same exponential distribution.
# Therefore, the expected time for your friend is the same as for any other customer.
expected_time_for_friend = expected_checkout_time_for_any_customer

print("The expected checkout time for the previous customer was given as {} minutes.".format(expected_checkout_time_for_any_customer))
print("Because the checkout process is the same for everyone, the expected time for your friend is also {} minutes.".format(expected_time_for_friend))