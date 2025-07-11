import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_covalency():
    """
    Analyzes the relative covalency of CeF6(2-) and CeCl6(2-)
    based on the provided information.
    """
    # Step 1: Define the core chemical principle.
    # Covalency is the degree of electron sharing between atoms in a chemical bond.
    # This sharing is enabled by the overlap of atomic orbitals.
    print("Step 1: Establishing the link between Orbital Overlap and Covalency.")
    print("In chemical bonding, the extent of orbital overlap between two atoms is directly proportional to the covalent character of the bond.")
    print("Greater orbital overlap allows for more effective sharing of electrons, which signifies stronger covalency.")
    print("-" * 60)

    # Step 2: Use the information given in the problem statement.
    # The problem explicitly states the relative overlap for the two compounds.
    print("Step 2: Applying the given experimental observation.")
    print("Given: The overlap of Ce(IV) 4f orbitals with F 2p orbitals in CeF6(2-) is GREATER than the overlap with Cl 3p orbitals in CeCl6(2-).")
    print("-" * 60)

    # Step 3: Draw the logical conclusion.
    # Based on the principle from Step 1 and the data from Step 2, a conclusion can be reached.
    print("Step 3: Drawing the conclusion.")
    print("Since CeF6(2-) has greater orbital overlap compared to CeCl6(2-), it follows directly that CeF6(2-) will display stronger covalency.")

# Execute the analysis
analyze_covalency()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the console
print(output)

# The final answer as a word.
final_answer = "stronger"