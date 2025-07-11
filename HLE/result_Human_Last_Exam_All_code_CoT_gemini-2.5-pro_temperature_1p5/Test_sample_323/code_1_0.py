import functools
import operator

def solve_nim_variant(scenarios):
    """
    Solves a series of modified Nim game scenarios.

    In this game, a "Zero-Move" is allowed, which passes the turn.
    Each of the n piles can have a Zero-Move applied to it at most once.

    The winning strategy depends on the number of piles (n) and the Nim-sum (bitwise XOR of pile sizes).
    - If n=1: The first player (Alice) can always clear the only pile and win.
    - If n > 1 and n is even: The game follows standard Nim rules. The first player wins if the Nim-sum is non-zero.
    - If n > 1 and n is odd: The game follows reversed Nim rules. The first player wins if the Nim-sum is zero.
    """
    final_result_string = ""
    
    for i, (n, a) in enumerate(scenarios, 1):
        print(f"Analyzing scenario ({i}): a={a}")
        print(f"  - Number of piles (n): {n}")
        
        winner = ''
        
        # Handle the special case where n=1
        if n == 1:
            nim_sum = a[0]
            equation_str = f"{a[0]}"
            print(f"  - Nim-sum = {equation_str} = {nim_sum}")
            winner = 'Alice'
            print(f"  - Rule: With n=1, Alice can always clear the single pile and win immediately.")
        
        # Handle the general case for n > 1
        else:
            nim_sum = functools.reduce(operator.xor, a)
            equation_str = ' ^ '.join(map(str, a))
            print(f"  - Nim-sum = {equation_str} = {nim_sum}")
            
            is_n_even = (n % 2 == 0)
            
            if is_n_even:
                # If n is even, it's like standard Nim
                if nim_sum != 0:
                    winner = 'Alice'
                else:
                    winner = 'Bob'
                print(f"  - Rule: With n={n} (even), the outcome follows standard Nim. Winner: {winner}")
            else:
                # If n is odd, the winning condition is reversed
                if nim_sum == 0:
                    winner = 'Alice'
                else:
                    winner = 'Bob'
                print(f"  - Rule: With n={n} (odd), the winning condition is reversed. Winner: {winner}")

        winner_char = winner[0]
        final_result_string += winner_char
        print(f"  --> Conclusion for scenario ({i}): {winner} wins.\n")

    print("="*30)
    print("Final result string:")
    print(final_result_string)


if __name__ == '__main__':
    # List of scenarios to solve
    # Format: (n, [a_1, a_2, ..., a_n])
    problem_scenarios = [
        (2, [12, 12]),
        (3, [8, 15, 7]),
        (3, [7, 16, 8]),
        (4, [12, 19, 21, 10]),
        (5, [16, 25, 58, 22, 60])
    ]
    
    solve_nim_variant(problem_scenarios)
    print("\n<<<BABBB>>>")