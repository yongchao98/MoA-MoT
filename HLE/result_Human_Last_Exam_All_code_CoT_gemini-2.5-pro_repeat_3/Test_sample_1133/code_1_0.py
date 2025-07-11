import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# --- Main Logic ---

# Step 1: Define the known information from the problem description.
# The problem states that the 10 minutes figure represents the expected value.
expected_checkout_time = 10  # in minutes

# Step 2: Explain the reasoning based on probability theory.
print("This problem relies on understanding the concept of 'Expected Value' for independent events.")
print(f"We are explicitly told that the expected checkout time for a customer is {expected_checkout_time} minutes.")
print("\nEach customer's checkout is an independent process. The duration of one customer's checkout does not influence the next.")

# Step 3: Apply the reasoning to the friend's situation.
# Since the friend's checkout is a new, independent instance of the same process,
# its expected duration is the same as the general expected duration.
expected_time_for_friend = expected_checkout_time

print("\nWhen your friend arrives at the cashier, a new, independent checkout process begins.")
print("The expected duration for this new process is the same as the average for any customer.")

# Step 4: Formulate and print the final answer equation.
print("\nTherefore, the expected time for your friend is calculated as follows:")
print(f"Expected Time (Friend) = Expected Time (Any Customer) = {expected_time_for_friend} minutes")

# --- End of Logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the content of the buffer
output = string_buffer.getvalue()

# Print the captured output
print(output)

# Final answer selection
final_answer = 'C'
print(f"\n<<<C>>>")