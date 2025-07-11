# The problem states that the expected value for a customer's checkout process is 10 minutes.
# We can assign this to a variable.
expected_value_of_checkout = 10

# Each customer's checkout is an independent event. The time taken by one customer
# does not influence the time taken by the next.
# Your friend is starting a new checkout process. Therefore, the expected time
# for your friend's checkout is the same as the given expected value for any customer.
friend_expected_time = expected_value_of_checkout

# The final equation is simply stating that the friend's expected time is equal
# to the known expected value of the process.
print(f"The expected time for any customer's checkout is given as {expected_value_of_checkout} minutes.")
print("Since the friend's checkout is a new, independent event from the same process,")
print(f"the expected time until the friend is done is also {friend_expected_time} minutes.")
print("\nFinal Equation:")
print(f"Friend's Expected Time = {friend_expected_time}")