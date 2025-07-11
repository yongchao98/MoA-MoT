import functools
import operator

def solve_modified_nim():
    """
    Solves the modified Nim game for a list of scenarios.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    final_result_string = ""

    print("Analyzing Nim game scenarios with Zero-Moves:")
    print("=" * 40)

    for i, piles in enumerate(scenarios, 1):
        n = len(piles)
        nim_sum = functools.reduce(operator.xor, piles)

        # Logic: Alice wins if (nim_sum == 0 and n is odd) or (nim_sum != 0 and n is even)
        # This is equivalent to (nim_sum == 0) != (n is even)
        is_nim_sum_zero = (nim_sum == 0)
        is_n_even = (n % 2 == 0)

        if (is_nim_sum_zero and not is_n_even) or (not is_nim_sum_zero and is_n_even):
            winner = "Alice"
            final_result_string += "A"
        else:
            winner = "Bob"
            final_result_string += "B"

        # Output the detailed analysis for the current scenario
        n_parity = "even" if is_n_even else "odd"
        equation = " ^ ".join(map(str, piles)) + f" = {nim_sum}"

        print(f"Scenario {i}: a = {piles}")
        print(f"  - Number of piles (n): {n} ({n_parity})")
        print(f"  - Nim-sum equation (S): {equation}")
        print(f"  - Verdict: {winner} wins.")
        print("-" * 20)

    print(f"Final combined result string: {final_result_string}")


solve_modified_nim()