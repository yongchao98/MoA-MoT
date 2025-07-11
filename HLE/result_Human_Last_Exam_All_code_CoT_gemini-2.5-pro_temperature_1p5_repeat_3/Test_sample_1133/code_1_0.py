# Step 1: Define the given information from the problem.
# The problem states that the checkout time of a previous customer was 10 minutes,
# and this represents the expected value (average) of the checkout process.
# This value is the parameter for our probability model.
expected_checkout_time = 10

# Step 2: Identify the process for your friend.
# Your friend is now starting their own checkout process. This is a new and
# independent instance of the same process. It follows the same rules and
# has the same average time.
# The question asks for the expected time until your friend is done.

# Step 3: State the conclusion.
# Because your friend's checkout follows the same probability distribution,
# their expected time is the same as the given average time.
expected_time_for_friend = expected_checkout_time

print("This problem is an application of the exponential distribution.")
print("The expected (average) checkout time for any customer is given.")
print(f"Expected_Time_Average = {expected_checkout_time} minutes")
print("\nYour friend's checkout is a new event that follows the same distribution.")
print("Therefore, their expected checkout time is the same as the average.")
print(f"Expected_Time_Friend = Expected_Time_Average")
print(f"\nExpected time for your friend: {expected_time_for_friend} minutes")
