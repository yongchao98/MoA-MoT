import sys
import io

# This problem is a classic example of applying probability concepts.
# Let's break it down using the information provided.

# 1. The problem describes a process with a constant probability of completion per unit of time.
#    This implies the checkout duration follows an exponential distribution.

# 2. We are explicitly given the expected value (the mean) of this distribution.
expected_time_for_any_customer = 10  # minutes

# 3. Your friend is a new customer. Their checkout process is a new, independent event
#    that follows the exact same exponential distribution as any other customer.

# 4. The question asks for the expected time for your friend's checkout.
#    Since it's the same process, the expected value is the same.
expected_time_for_friend = expected_time_for_any_customer

# The information about you waiting for the previous customer is extra detail
# to set the scene but doesn't change the expected outcome for your friend's
# separate and new checkout process.

print("This problem tests the understanding of independent events in a statistical process.")
print(f"The given expected checkout time for any customer is: {expected_time_for_any_customer} minutes.")
print("Each customer's checkout is an independent event drawn from the same distribution.")
print("Therefore, your friend's expected checkout time is the same as the given expected time.")

print("\nFinal Equation:")
# The "equation" here is a statement of identity, as the logic dictates the values are the same.
print(f"Expected_Time_Friend = Expected_Time_Any_Customer")
print(f"{expected_time_for_friend} minutes = {expected_time_for_any_customer} minutes")

print(f"\nThe expected time until your friend is done is {expected_time_for_friend} minutes.")
