import sys
import io

# Keep the original standard output
original_stdout = sys.stdout
# Redirect standard output to a new stream
sys.stdout = captured_output = io.StringIO()

# Design a molecule that fulfills the given criteria.
# The molecular formula is determined to be C12H24N2O3 based on the molecular weight and valence electron count.
# This formula allows for 3 ether oxygens and 2 tertiary amines, fitting within the 5 total heteroatom limit.
# The contradiction in the prompt (requesting 5 ether oxygens with a formula that only supports 3) is resolved
# by assuming a typo and proceeding with 3 ether oxygens.
# The structure consists of an N,N'-dimethylpiperazine ring connected to a 1,3,5-trioxane ring via a propyl linker.

# SMILES representation of 2-(3-(N,N'-dimethylpiperazin-2-yl)propyl)-1,3,5-trioxane
smiles_representation = "CN1CCN(C)C(C1)CCCC2OCOCO2"

# Print the final SMILES string
print(smiles_representation)

# Restore the original standard output
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the output to the user and also add the final answer tag
print(output, end="")
print(f"<<<{output.strip()}>>>")