import functools
import operator

def solve_nim_variant():
    """
    Solves the modified Nim game for a given set of scenarios.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]

    final_result_chars = []

    print("Analyzing Nim game scenarios with Zero-Moves:\n")

    for i, s in enumerate(scenarios):
        n = s['n']
        a = s['a']
        
        # Calculate the Nim-sum
        nim_sum = functools.reduce(operator.xor, a)
        
        # Format the equation string as requested
        nim_sum_equation = " XOR ".join(map(str, a))

        print(f"--- Scenario ({i+1}) ---")
        print(f"Piles: a = {a}")
        print(f"Number of piles: n = {n}")
        print(f"Nim-sum calculation: {nim_sum_equation} = {nim_sum}")

        winner = ''
        explanation = ''

        if n % 2 == 0:  # n is even
            if nim_sum != 0:
                winner = 'A'
                explanation = "n is even and Nim-sum is non-zero. Alice wins."
            else:
                winner = 'B'
                explanation = "n is even and Nim-sum is zero. Bob wins."
        else:  # n is odd
            if nim_sum == 0:
                winner = 'A'
                explanation = "n is odd and Nim-sum is zero. Alice wins."
            else:
                winner = 'B'
                explanation = "n is odd and Nim-sum is non-zero. Bob wins."
        
        print(f"Result: {explanation}\n")
        final_result_chars.append(winner)

    final_answer = "".join(final_result_chars)
    print("=========================")
    print(f"Final combined result string: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_nim_variant()