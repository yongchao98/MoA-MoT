import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The chosen statements are II and III.
# Statement II provides a practical first-line strategy for cessation (NRT).
# Statement III provides the core reason why an adolescent must quit (unknown risks).
# Together, they form the most appropriate counseling advice.

statement_II = 2
statement_III = 3

print(f"The best options to consider are {statement_II} and {statement_III}.")
print("This combination corresponds to answer choice J.")

# Get the captured output
response_text = captured_output.getvalue()
# Restore original stdout
sys.stdout = original_stdout

# Final print to the user
print(response_text, end="")
print("<<<J>>>")