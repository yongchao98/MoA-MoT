# The problem states that the expected checkout time for a customer is 10 minutes.
# This value represents the mean of the underlying probability distribution (an exponential distribution).
expected_time_for_a_customer = 10

# Each customer's checkout is an independent event drawn from this same distribution.
# Therefore, the expected checkout time for your friend is identical to the expected time for any other customer.
expected_time_for_friend = expected_time_for_a_customer

# The final equation is simply that the expected time for the friend is equal to the given expected time for any customer.
# Equation: E(friend's time) = 10
# The number in this final equation is 10.
print("The final equation is: Expected Time = 10 minutes.")
print("The number in the final equation is:")
print(expected_time_for_friend)