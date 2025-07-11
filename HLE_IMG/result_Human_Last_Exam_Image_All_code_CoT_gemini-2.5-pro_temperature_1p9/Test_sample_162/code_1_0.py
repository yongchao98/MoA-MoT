import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

def identify_reactant():
  """
  This function identifies the missing reactant in the given chemical synthesis.
  The synthesis involves a Wittig reaction followed by a Michael-initiated ring-closure (MIRC) reaction.
  The MIRC part transforms an alpha,beta-unsaturated ketone into a cyclohexane-1,3-dione.
  This specific transformation is achieved by using a malonic ester as the Michael donor.
  The reaction sequence involves:
  1. Michael addition of the malonate enolate.
  2. Intramolecular Claisen condensation to form the 6-membered ring.
  3. Saponification of the remaining ester group.
  4. Decarboxylation of the resulting beta-keto acid.
  This analysis points to diethyl malonate as the required reactant.
  """
  reactant_name = "diethyl malonate"
  print(reactant_name)

identify_reactant()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the captured output to the console
print(output)