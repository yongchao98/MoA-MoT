import functools
import operator

def solve_nim_scenarios():
    """
    Solves a series of modified Nim game scenarios.

    In this game, a "Zero-Move" allows a player to pass their turn,
    using up the Zero-Move for a specific pile. The winner is determined
    by the Nim-sum (XOR sum of pile sizes) and the number of piles (n).

    The logic is as follows:
    1. If the Nim-sum is non-zero, the first player (Alice) has a winning
       strategy from standard Nim and wins.
    2. If the Nim-sum is zero, the game becomes a "passing contest" using
       the n available Zero-Moves. The winner is determined by who gets
       the last pass. If n is odd, Alice wins. If n is even, Bob wins.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]},
    ]

    final_result_string = ""
    print("--- Nim Game Analysis ---")

    for i, case in enumerate(scenarios):
        n = case['n']
        a = case['a']
        
        # Calculate the Nim-sum of the piles
        nim_sum = functools.reduce(operator.xor, a)
        
        equation_str = " ^ ".join(map(str, a))
        print(f"\nScenario ({i+1}): n={n}, a={a}")
        print(f"Nim-sum equation: {equation_str} = {nim_sum}")
        
        winner = ''
        if nim_sum != 0:
            winner = 'A'
            print(f"Result: Nim-sum is non-zero. Alice wins.")
        else:  # nim_sum == 0
            if n % 2 != 0:
                winner = 'A'
                print(f"Result: Nim-sum is 0 and n={n} is odd. Alice wins.")
            else:  # n is even
                winner = 'B'
                print(f"Result: Nim-sum is 0 and n={n} is even. Bob wins.")
        
        final_result_string += winner

    print("\n--- Final Answer ---")
    print("The concatenated string of winners (A for Alice, B for Bob) is:")
    print(final_result_string)

solve_nim_scenarios()