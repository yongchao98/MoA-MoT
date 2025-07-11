# The problem describes a checkout process where the completion probability is constant
# over time. This implies the checkout duration follows an exponential distribution.
# A key feature of this distribution is the "memoryless property".

# This property means that the expected future waiting time is independent of the time
# that has already passed. Each customer's checkout is an independent event.

# We are given a crucial piece of information: the expected checkout time for a customer.
# The problem states we can assume the 10 minutes is this expected value.
expected_checkout_time = 10

# Your friend is starting their checkout process now. Since their checkout is a new
# and independent event drawn from the same distribution, their expected time
# is the same as the expected time for any other customer.

expected_time_for_friend = expected_checkout_time

# The final equation is straightforward:
# E(Friend's Time) = E(Any Customer's Time)
# where E() denotes the expected value.

print("Let E_customer be the expected checkout time for any customer.")
print("From the problem, we are given the value for this expectation.")
print(f"E_customer = {expected_checkout_time} minutes")
print("\nLet E_friend be the expected checkout time for your friend.")
print("Since the friend's checkout is a new, independent process, it has the same expected value.")
print("E_friend = E_customer")
print(f"E_friend = {expected_time_for_friend}")
print(f"\nTherefore, the expected time until your friend is done is {expected_time_for_friend} minutes.")