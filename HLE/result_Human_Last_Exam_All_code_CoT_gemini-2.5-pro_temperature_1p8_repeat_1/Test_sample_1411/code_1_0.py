import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# This is a knowledge-based question. Research shows the watch is a
# Rolex Submariner, Ref. 5513, that was gifted by Steve McQueen.
# It was uniquely customized by George "Stan" Shourds, who added a
# date function and his own name to the dial.

# The text located directly above the date window.
text_on_watch = "shourds"

# The user requested the answer in all lower case.
final_answer = text_on_watch.lower()

# Print the final answer to the captured output.
# The user's prompt has an unusual instruction: "Remember in the final code
# you still need to output each number in the final equation!".
# Since there is no equation or number, this instruction is not applicable,
# and we will print the requested text answer.
print(final_answer)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the captured output to the actual console
print(output)