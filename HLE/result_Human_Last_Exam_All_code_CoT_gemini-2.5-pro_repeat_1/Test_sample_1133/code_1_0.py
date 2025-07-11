# The problem states that the expected checkout time for a customer is 10 minutes.
# This value defines the average rate of the checkout process.
expected_time_per_customer = 10 # in minutes

# Let's denote the expected time for your friend as E[Friend].
# The checkout process is the same for every customer. Each customer's checkout time
# is an independent draw from the same probability distribution.

# Therefore, the expected checkout time for your friend is the same as the
# established expected time for any customer.
expected_time_for_friend = expected_time_per_customer

# The final equation is simply stating this equality.
print(f"The expected time for any customer is {expected_time_per_customer} minutes.")
print("Since the friend's checkout is an independent event from the same process,")
print(f"the expected time for the friend is also {expected_time_for_friend} minutes.")
print("\nFinal Equation:")
print(f"Expected_Time(Friend) = Expected_Time(Any_Customer)")
print(f"Expected_Time(Friend) = {expected_time_for_friend}")
