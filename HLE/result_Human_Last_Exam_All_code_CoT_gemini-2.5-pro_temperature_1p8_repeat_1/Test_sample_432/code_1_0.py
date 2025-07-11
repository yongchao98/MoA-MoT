def solve_cardinality_puzzle():
    """
    Identifies which of a given list of sets have the same cardinality as the interval [0, 1].

    The cardinality of [0, 1] is the cardinality of the continuum, c, where c = 2^aleph_0.
    This script evaluates the cardinality of each set based on principles of set theory
    and prints the sorted string of letters corresponding to the sets with cardinality c.
    """
    
    # A dictionary to hold the cardinality of each set, using 'c' for continuum,
    # 'aleph_0' for countable infinity, and '>c' for cardinalities greater than continuum.
    cardinalities = {
        'A': 'c',   # (0, 1) has a bijection to R, so its cardinality is c.
        'B': 'aleph_0', # The natural numbers N are countably infinite.
        'C': 'aleph_0', # The rational numbers Q are countably infinite.
        'D': 'c',   # The real numbers R have cardinality c by definition.
        'E': 'c',   # The irrational numbers R \ Q. Since R = Q U (R \ Q), c = aleph_0 + |R \ Q|, so |R \ Q| = c.
        'F': 'c',   # The complex numbers C are equivalent to R^2. |R^2| = c^2 = c.
        'G': 'c',   # The quaternions H are equivalent to R^4. |R^4| = c^4 = c.
        'H': 'c',   # The set {x: c'(x) = 0} for the Cantor function contains open intervals, thus has cardinality c.
        'I': 'aleph_0', # The set of all finite strings over a finite/countable alphabet is countably infinite.
        'J': 'c',   # R^N, the set of all real-valued sequences. |R^N| = c^aleph_0 = c.
        'K': 'c',   # Z^N, the set of all integer sequences. |Z^N| = aleph_0^aleph_0 = c.
        'L': 'aleph_0', # N x N is a Cartesian product of countable sets, hence countable.
        'M': 'c',   # R x R is equivalent to R^2, with cardinality c^2 = c.
        'N': 'c',   # 2^N, the power set of naturals. |2^N| = 2^aleph_0 = c.
        'O': 'c',   # 2^Q, the power set of rationals. |Q|=aleph_0, so |2^Q|=2^aleph_0 = c.
        'P': '>c',  # 2^C, the power set of complex numbers. |C|=c, so |2^C| = 2^c, which is > c.
        'Q': '>c'   # 2^(2^Z). |Z|=aleph_0, so cardinality is 2^(2^aleph_0) = 2^c, which is > c.
    }
    
    # Find all options with cardinality 'c'.
    correct_options = []
    for option, card in cardinalities.items():
        if card == 'c':
            correct_options.append(option)
            
    # Sort the results alphabetically.
    correct_options.sort()
    
    # Join into a single string and print.
    final_answer = "".join(correct_options)
    
    print(final_answer)

if __name__ == '__main__':
    solve_cardinality_puzzle()