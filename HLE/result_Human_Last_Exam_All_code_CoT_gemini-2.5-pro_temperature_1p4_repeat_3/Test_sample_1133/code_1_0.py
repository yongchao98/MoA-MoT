# The problem describes a checkout process where the probability of completion
# is constant over time. This scenario is best modeled by the exponential distribution,
# which has a key property: memorylessness.

# The memoryless property means that the expected remaining time for an event
# is independent of how much time has already passed.

# We are told that the expected value of a customer's total checkout time is 10 minutes.
expected_time_per_customer = 10

# Your friend is starting their checkout process from the beginning.
# Therefore, their expected checkout time is the same as the average checkout time
# for any customer, which is given as 10 minutes.
expected_time_for_friend = expected_time_per_customer

# The final calculation is a direct application of this principle.
print("The final equation is:")
# We show that the expected time for the friend is equal to the known expected time.
print(f"{expected_time_for_friend} = {expected_time_per_customer}")

print("\nTherefore, the expected time until your friend is done and ready to go is:")
print(f"{expected_time_for_friend} minutes.")