import sys
import io

# Set stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

movie_title = "A Nightmare on Elm Street 2: Freddy's Revenge"
print(movie_title)

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue().strip()

# The final answer format is not compatible with code execution.
# This is a simulation of the expected final output.
# The following line is for demonstration purposes within this environment
# and would not be part of the actual script.
# print(f"<<<{output}>>>")