import math

# The problem states that the average, or expected, checkout time for a customer is 10 minutes.
# This is the mean of the underlying probability distribution (an exponential distribution).
expected_checkout_time_per_customer = 10

# The checkout process is identical and independent for every customer.
# Your friend is a new customer starting their checkout process from the beginning.
# Therefore, the expected time for your friend's checkout is the same as the
# average time for any other customer.

expected_time_for_friend = expected_checkout_time_per_customer

# The problem is a demonstration of the memoryless property. The time it took
# for a previous customer simply establishes the average for the process.
# It doesn't influence the duration of the next independent process.

print("Equation of Expected Time:")
print(f"Expected_Time_for_Friend = Expected_Checkout_Time_Per_Customer")
print(f"Expected_Time_for_Friend = {expected_checkout_time_per_customer}")
print(f"\nThe expected time until your friend is done is {expected_time_for_friend} minutes.")
