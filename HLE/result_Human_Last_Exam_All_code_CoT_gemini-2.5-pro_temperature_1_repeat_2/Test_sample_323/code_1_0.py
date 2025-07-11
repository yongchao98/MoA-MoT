import functools
import operator

def solve_nim_scenarios():
    """
    Solves a series of modified Nim game scenarios.

    The winning condition depends on the parity of the number of piles (n)
    and whether the Nim-sum (XOR sum) of the piles is zero.

    - If n is even:
        - Nim-sum != 0 -> Alice wins
        - Nim-sum == 0 -> Bob wins
    - If n is odd:
        - Nim-sum == 0 -> Alice wins
        - Nim-sum != 0 -> Bob wins
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    final_answer = ""

    for i, piles in enumerate(scenarios):
        n = len(piles)
        
        # Calculate the Nim-sum using functools.reduce and operator.xor
        nim_sum = functools.reduce(operator.xor, piles)

        # Build the equation string for a clear output
        equation_str = " XOR ".join(map(str, piles)) + f" = {nim_sum}"
        
        is_n_even = (n % 2 == 0)
        is_nim_sum_zero = (nim_sum == 0)
        
        winner = ''
        explanation = ''
        
        # Determine the winner based on the derived logic
        if is_n_even:
            if is_nim_sum_zero:
                winner = 'B'
                explanation = f"n={n} (even) and Nim-sum is 0. Bob wins."
            else:
                winner = 'A'
                explanation = f"n={n} (even) and Nim-sum is non-zero. Alice wins."
        else:  # n is odd
            if is_nim_sum_zero:
                winner = 'A'
                explanation = f"n={n} (odd) and Nim-sum is 0. Alice wins."
            else:
                winner = 'B'
                explanation = f"n={n} (odd) and Nim-sum is non-zero. Bob wins."
        
        print(f"Scenario ({i+1}): a={piles}")
        print(f"Calculation: {equation_str}")
        print(f"Result: {explanation}")
        print("-" * 30)
        
        final_answer += winner
        
    print(f"The final combined result is: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_nim_scenarios()