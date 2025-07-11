import math

# The problem states that the expected checkout time for a customer is 10 minutes.
# This is the mean of the underlying exponential distribution.
# Let's represent this given expected value.
expected_time_for_any_customer = 10

# The checkout process is the same for every customer.
# Therefore, the expected checkout time for your friend is the same as for any other customer.
expected_time_for_friend = expected_time_for_any_customer

print("The checkout process is modeled by an exponential distribution due to the constant probability of completion per second.")
print("A key property of this distribution is that each event (each customer's checkout) is independent.")
print(f"The expected checkout time for any customer is given as: {expected_time_for_any_customer} minutes.")
print(f"Since your friend's checkout is a new, independent event from the same distribution, their expected time is also: {expected_time_for_friend} minutes.")
print(f"Final Equation: E(friend's time) = {expected_time_for_friend}")
