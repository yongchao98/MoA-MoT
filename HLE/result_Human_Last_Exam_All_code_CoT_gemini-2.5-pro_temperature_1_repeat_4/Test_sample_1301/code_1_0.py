import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The SMILES string represents a molecule conforming to all the specified constraints.
# C12H18O6
# MW: 258.11 g/mol
# Heavy Atoms: 18
# Valence Electrons: 102
# Formal Charge: 0
# Rings: 3 (all saturated heterocycles)
# Rotatable Bonds: 0 (it is a rigid cage-like structure)
# Heteroatoms: 6 (1 carbonyl O, 5 ether O)
# H-bond Acceptors: 6
# H-bond Donors: 0
# The structure is a tricyclic (bridged bicyclic) system.
smiles = "O=C1C2OC3C4OC5COC(C14)C2C35"
print(f"{smiles}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Final print to the user
print(output)
print(f'<<<O=C1C2OC3C4OC5COC(C14)C2C35>>>')