import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The user wants to identify the name of the reactant.
# Based on the chemical transformation, which is a Michael addition followed by a Dieckmann cyclization,
# hydrolysis, and decarboxylation to form a 5-substituted cyclohexane-1,3-dione from an alpha,beta-unsaturated ketone,
# the required reactant is a malonic ester. The most common one used is diethyl malonate.

reactant_name = "Diethyl malonate"
print(reactant_name)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)