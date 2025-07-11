# The problem describes a scenario where the time to complete a task (checkout)
# follows a specific probability distribution.

# Step 1: Identify the given information.
# We are told that the expected checkout time for a customer is 10 minutes.
expected_checkout_time = 10

# Step 2: Understand the process for different customers.
# The problem states that the process "is the same for everyone".
# This means that the checkout time for each customer is an independent random event
# drawn from the same probability distribution.
# Let T_customer be the checkout time for the first customer.
# Let T_friend be the checkout time for your friend.
# We are given E[T_customer] = 10 minutes.
# We need to find E[T_friend].

# Step 3: Apply the principle of identical distribution.
# Since both checkout times are drawn from the same distribution, their
# expected values (averages) must be the same.
expected_time_for_friend = expected_checkout_time

# Step 4: Print the final result and the reasoning.
print("The expected checkout time for one customer is given as 10 minutes.")
print("The problem states the process is the same for all customers.")
print("Therefore, your friend's checkout process is an independent event with the same expected duration.")
print("The expected time until your friend is done is:")
print(f"{expected_time_for_friend}")