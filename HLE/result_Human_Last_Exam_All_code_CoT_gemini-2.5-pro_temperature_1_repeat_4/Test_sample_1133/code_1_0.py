# The problem states that the checkout time for a customer follows a specific probability distribution.
# The key phrase "probability that a person completes the checkout process is constant per second"
# implies that the duration of the checkout follows an exponential distribution.

# We are given that the expected value (the average time) of this process is 10 minutes,
# based on the information from a previous customer.
expected_time_for_any_customer = 10  # in minutes

# Your friend's checkout is a new, independent event that follows the same distribution.
# The expected value for any single event from this distribution is the same.
expected_time_for_friend = expected_time_for_any_customer

# The final "equation" is simply stating that the expected time for your friend
# is equal to the known expected time of the process.
print(f"The given expected checkout time for a customer is: {expected_time_for_any_customer} minutes.")
print(f"Since your friend's checkout is a new and independent event from the same process,")
print(f"the expected time until your friend is done is also: {expected_time_for_friend} minutes.")
print(f"The final logical step is: Expected_Time_Friend ({expected_time_for_friend}) = Expected_Time_Process ({expected_time_for_any_customer})")
