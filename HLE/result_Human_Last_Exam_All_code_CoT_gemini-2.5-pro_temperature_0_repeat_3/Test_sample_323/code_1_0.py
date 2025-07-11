import functools
import operator

def solve_nim_scenarios():
    """
    Solves the modified Nim game for a list of predefined scenarios.
    For each scenario, it calculates the Nim-sum and applies the winning logic
    based on the Nim-sum and the number of piles (n).
    It prints the analysis for each case and the final result string.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    final_result_string = ""

    for i, piles in enumerate(scenarios):
        n = len(piles)
        # Calculate the Nim-sum using XOR operation on all pile sizes
        nim_sum = functools.reduce(operator.xor, piles)

        # Format the equation string as requested
        equation_str = " XOR ".join(map(str, piles))
        
        print(f"Scenario ({i+1}) with piles {piles}:")
        print(f"Equation: {equation_str} = {nim_sum}")

        winner = ''
        # Case 1: Nim-sum is non-zero (standard Nim winning position)
        if nim_sum != 0:
            # If n is even, Alice wins. If n is odd, Bob wins.
            if n % 2 == 0:
                winner = 'A'
                print(f"Result: Nim-sum is non-zero and n={n} is even. Alice wins.")
            else:
                winner = 'B'
                print(f"Result: Nim-sum is non-zero and n={n} is odd. Bob wins.")
        # Case 2: Nim-sum is zero (standard Nim losing position)
        else:
            # If n is odd, Alice wins. If n is even, Bob wins.
            if n % 2 != 0:
                winner = 'A'
                print(f"Result: Nim-sum is zero and n={n} is odd. Alice wins.")
            else:
                winner = 'B'
                print(f"Result: Nim-sum is zero and n={n} is even. Bob wins.")
        
        final_result_string += winner
        print("-" * 30)

    print(f"Final answer string: {final_result_string}")

# Execute the function to solve the scenarios
solve_nim_scenarios()