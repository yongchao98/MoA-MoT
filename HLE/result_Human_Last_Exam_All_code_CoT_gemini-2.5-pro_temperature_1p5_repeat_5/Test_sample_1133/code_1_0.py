# The problem describes a process with a constant probability of completion over time,
# which is characteristic of an exponential distribution. A key feature of this
# distribution is that it is "memoryless," meaning past events don't influence
# future outcomes for independent trials.

# We are explicitly told that the expected checkout time for a customer is 10 minutes.
# This establishes the mean of the distribution for this process.
expected_checkout_time_minutes = 10

# Your friend is starting their checkout process. Since the process is the same
# for every customer, your friend's checkout is a new, independent event
# drawn from the same distribution.
# Therefore, the expected time for your friend is the same as the established
# mean of the distribution.
expected_time_for_friend = expected_checkout_time_minutes

# The final equation is:
print(f"E[Friend's checkout time] = {expected_time_for_friend} minutes")