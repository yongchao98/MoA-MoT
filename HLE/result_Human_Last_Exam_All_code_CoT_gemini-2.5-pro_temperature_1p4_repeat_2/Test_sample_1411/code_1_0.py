import sys
import io

# Set stdout to a string buffer to capture output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The watch is the Heuer Monaco, Ref. 1133, worn by Steve McQueen in 'Le Mans'.
# The date window is at the 6 o'clock position.
# At the 3 o'clock position, the words 'automatic' and 'chronograph' are printed one above the other.
# The word 'chronograph' is the lower of the two and is the text positioned directly above the date window.
# The user requested the answer in all lower case.
the_watch_text = "chronograph"

print(the_watch_text)

# Restore stdout and get the captured content
sys.stdout = old_stdout
output = captured_output.getvalue().strip()

# Final answer structure
# No need to print this part in the final execution. This is for structure.
# print("<<<{}>>>".format(output))

# Final execution print
print(output)