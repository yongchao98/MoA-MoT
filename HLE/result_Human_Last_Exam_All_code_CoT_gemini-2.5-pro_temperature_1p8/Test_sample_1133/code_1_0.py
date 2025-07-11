# The problem describes a process where the probability of completion per unit of time is constant.
# This implies the checkout time follows an exponential distribution, which is memoryless.

# We are given the expected value (mean) of this distribution from a previous customer's experience.
# The prompt states to assume the "10 minutes" is the expected value.
mean_checkout_time = 10

# The friend's checkout process is governed by the same distribution, as it is "the same for everyone."
# The question asks for the expected time for the friend's checkout process, starting from when they
# reach the cashier.
# For an identical process, the expected value remains the same.
expected_time_for_friend = mean_checkout_time

# The final equation is simply stating that the friend's expected time is equal to the given mean time.
print(f"Given mean checkout time = {mean_checkout_time} minutes")
print(f"Expected time for friend = {expected_time_for_friend} minutes")
