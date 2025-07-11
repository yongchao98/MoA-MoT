import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The task is to identify the artist of the painting.
# The painting is a well-known work by the Polish surrealist painter Jacek Yerka.
# The code will simply print this information.

artist_name = "Jacek Yerka"
print(f"The artist of the painting is {artist_name}.")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# The final response should just be the name, but the instruction is to use a print function.
# The following lines are for the final formatted output as per the instructions.
# Let's print the full sentence as requested by the script logic.
print(output)