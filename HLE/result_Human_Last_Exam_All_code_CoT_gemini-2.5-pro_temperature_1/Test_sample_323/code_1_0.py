import functools
import operator

def solve_nim_with_zero_moves():
    """
    Solves a series of modified Nim game scenarios and prints the outcome.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    final_result = []

    for i, scenario in enumerate(scenarios):
        n = scenario['n']
        a = scenario['a']
        
        # Calculate the Nim-sum of the piles
        nim_sum = functools.reduce(operator.xor, a)
        
        # Build the equation string for the Nim-sum calculation
        equation_str = " XOR ".join(map(str, a)) + f" = {nim_sum}"
        
        n_is_even = (n % 2 == 0)
        
        winner = ''
        winner_initial = ''

        if n_is_even:
            # If n is even, standard Nim rules apply for the winner
            if nim_sum != 0:
                winner = "Alice"
                winner_initial = "A"
            else:
                winner = "Bob"
                winner_initial = "B"
        else:  # n is odd
            # If n is odd, the winner is the opposite of standard Nim
            if nim_sum == 0:
                winner = "Alice"
                winner_initial = "A"
            else:
                winner = "Bob"
                winner_initial = "B"
        
        print(f"({i+1}) n={n}, a={a}")
        print(f"    Nim-sum: {equation_str}")
        print(f"    n is {'even' if n_is_even else 'odd'}. According to the rules, {winner} wins.")
        
        final_result.append(winner_initial)
        
    print("\nFinal answer string:")
    print("".join(final_result))

solve_nim_with_zero_moves()