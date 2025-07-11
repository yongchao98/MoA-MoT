import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a new stream to capture the output
captured_output = io.StringIO()
# Redirect stdout to the new stream
sys.stdout = captured_output

# --- Start of the logic ---

# The problem states that the expected checkout time for a customer is 10 minutes.
# This is the mean of the underlying probability distribution.
expected_value = 10

# The key insight is that a process with a constant probability of completion per unit of time
# follows an exponential distribution. This distribution has a "memoryless" property.

# The memoryless property means that the expected waiting time for an event to occur is
# constant, regardless of how much time has already elapsed.

# Your friend is starting a new, independent checkout process. Therefore, their expected
# checkout time is simply the mean of the distribution, which was given.

# The equation for the expected time for your friend is:
# Expected_Time_Friend = Expected_Time_Any_Customer
expected_time_friend = expected_value

print("The scenario described is modeled by an exponential distribution, which is 'memoryless'.")
print("This means the expected waiting time for a new event is always the same, regardless of past events.")
print("The expected checkout time for any customer is given as 10 minutes.")
print("\nTherefore, the final equation for your friend's expected checkout time is:")
print(f"Expected_Time = {expected_time_friend}")
print("\nThe expected time until your friend is done is 10 minutes.")

# --- End of the logic ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output = captured_output.getvalue()

# Print the captured output
print(output)