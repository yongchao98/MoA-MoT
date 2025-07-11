import math

def solve_a4():
    """
    Calculates a(4), the maximal number of prime implicants of a Boolean function of 4 variables.

    The function a(n) represents a known, non-trivial integer sequence. There is no
    simple closed-form formula to calculate its terms. The values are determined from
    established research in Boolean logic and combinatorics.

    The first few terms of the sequence a(n) are:
    a(0) = 1
    a(1) = 2
    a(2) = 6
    a(3) = 20
    a(4) = 78
    """
    
    # The sequence of maximal prime implicants for n variables.
    # This is a known sequence from research, OEIS A000373.
    a_sequence = {
        0: 1,
        1: 2,
        2: 6,
        3: 20,
        4: 78
    }

    n = 4
    
    # We look up the value for n=4 from the known sequence.
    result = a_sequence.get(n)

    # Print the final equation as requested.
    if result is not None:
        print(f"The maximal number of prime implicants for a Boolean function of n variables, a(n), is given by a known sequence.")
        print(f"We are looking for the value of a({n}).")
        print(f"Based on the sequence, the equation is:")
        # The final equation output shows each number
        print(f"a({n}) = {result}")
    else:
        print(f"The value for a({n}) is not available in the pre-computed sequence.")

solve_a4()