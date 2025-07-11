import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a buffer
sys.stdout = captured_output = io.StringIO()

# The problem states that the expected (average) time for a customer's checkout is 10 minutes.
# This is the key parameter we are given.
expected_checkout_time = 10

# Your friend starts their checkout process. This is a new, independent instance of the same process.
# Therefore, the expected time for your friend to finish is the same as the given average.
expected_time_for_friend = expected_checkout_time

# The final equation simply expresses that the friend's expected time is this given value.
print(f"Let E_friend be the expected time for the friend.")
print(f"Let E_avg be the average checkout time given in the problem.")
print(f"E_friend = E_avg")
print(f"Therefore, the expected time for the friend is: {expected_time_for_friend} minutes.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)

# The final answer is one of the choices.
# Based on the calculation, the answer is 10 minutes, which corresponds to option C.
final_answer = "C"
print(f"\n<<<C>>>")