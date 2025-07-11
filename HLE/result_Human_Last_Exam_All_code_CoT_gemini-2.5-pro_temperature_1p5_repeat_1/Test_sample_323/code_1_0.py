import functools
import operator

def solve_nim_variant():
    """
    Solves a series of modified Nim game scenarios.

    In this variant of Nim, a "Zero-Move" is allowed, which acts as a pass.
    Each of the 'n' piles can be the subject of a Zero-Move once.

    The winning strategy depends on two factors:
    1. The Nim-sum (bitwise XOR sum) of the pile sizes.
    2. The parity (even or odd) of the number of piles, 'n', which is the
       initial number of available Zero-Moves.

    The logic is as follows:
    - Let S be the Nim-sum of all piles.
    - Let n be the number of piles.

    - If S is non-zero (a standard winning position for the first player):
      The first player (Alice) should make a standard Nim move to make S=0.
      The second player (Bob) must now use a Zero-Move to avoid losing immediately.
      They will then trade Zero-Moves. The winner of this sub-game depends on who
      uses the last Zero-Move.
      - If n is even, Alice wins.
      - If n is odd, Bob wins.

    - If S is zero (a standard losing position for the first player):
      Alice's only good move is to use a Zero-Move.
      They will then trade Zero-Moves.
      - If n is odd, Alice uses the last Zero-Move and wins.
      - If n is even, Bob uses the last Zero-Move and wins.

    Summary of winning conditions for Alice (first player):
    - Alice wins if (S != 0 and n is even) OR (S == 0 and n is odd).
    - Otherwise, Bob wins.
    """
    scenarios = [
        {"n": 2, "a": [12, 12]},
        {"n": 3, "a": [8, 15, 7]},
        {"n": 3, "a": [7, 16, 8]},
        {"n": 4, "a": [12, 19, 21, 10]},
        {"n": 5, "a": [16, 25, 58, 22, 60]}
    ]

    result_string = ""

    for scenario in scenarios:
        n = scenario["n"]
        piles = scenario["a"]

        # Calculate the Nim-sum (XOR sum) of the piles
        nim_sum = functools.reduce(operator.xor, piles)

        # Determine if n is even
        n_is_even = (n % 2 == 0)
        
        # Determine if the Nim-sum is zero
        nim_sum_is_zero = (nim_sum == 0)

        # Alice wins if (S!=0 and n is even) or (S==0 and n is odd)
        # This is equivalent to (nim_sum_is_zero != n_is_even)
        if nim_sum_is_zero != n_is_even:
            result_string += "A"
        else:
            result_string += "B"

    print(result_string)

solve_nim_variant()
<<<BABBB>>>