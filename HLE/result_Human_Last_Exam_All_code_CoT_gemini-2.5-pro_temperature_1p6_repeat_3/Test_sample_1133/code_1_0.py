import math

# The problem states that the probability of completing the checkout is constant per second.
# This implies the time 'T' for the checkout process follows an exponential distribution.
# A key property of the exponential distribution is that it is 'memoryless'.

# We are given that the expected time for a customer's checkout is 10 minutes.
# For an exponential distribution, the expected value E[T] is the mean of the distribution.
expected_time = 10

# Your friend is starting a new, independent checkout process.
# Because the process is memoryless, the expected time for your friend is simply the
# mean of the distribution, regardless of what happened before.
expected_time_for_friend = expected_time

# The question asks for the expected time until your friend is done from the moment they arrive at the cashier.
# Final Equation: Expected_Time_Friend = Expected_Time_Any_Customer
print("Let E_friend be the expected time for the friend.")
print("Let E_customer be the given expected time for any customer.")
print("Due to the memoryless property of the underlying exponential distribution, the two are equal.")
print("\nFinal Equation:")
print(f"E_friend = E_customer")
print(f"E_friend = {expected_time_for_friend}")

print(f"\nThe expected time until your friend is done is {expected_time_for_friend} minutes.")