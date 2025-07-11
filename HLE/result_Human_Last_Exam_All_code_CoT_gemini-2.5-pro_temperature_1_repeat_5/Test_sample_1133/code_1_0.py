# The problem states that the expected checkout time for a customer is 10 minutes.
# This value is derived from observing a previous customer.
expected_time_observed = 10  # in minutes

# --- Reasoning ---
# 1. The problem describes a process with a constant probability of completion per unit of time.
#    This means the checkout duration follows an exponential distribution.
# 2. The expected value (mean) of this distribution is given as 10 minutes.
# 3. Each customer's checkout is an independent trial. The duration for one customer
#    does not influence the duration for the next.
# 4. Your friend is starting a new, independent checkout process. This process is drawn
#    from the same distribution as the previous customer's.

# Therefore, the expected time for your friend's checkout is the same as the
# given expected time for any customer.

expected_time_for_friend = expected_time_observed

# --- Final Equation ---
# The logic is a direct equality based on the independent nature of the events.
# Expected Time (Friend) = Expected Time (Any Customer)

print(f"The given expected checkout time for a customer is: {expected_time_observed} minutes.")
print("Since each customer's checkout is an independent event from the same process,")
print("the expected time for the friend is the same.")
print(f"Final Equation: Expected Time (Friend) = {expected_time_for_friend}")
print(f"The expected time until your friend is done is {expected_time_for_friend} minutes.")
