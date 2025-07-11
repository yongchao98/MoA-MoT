import functools
import operator

def solve_nim_variant():
    """
    Solves the modified Nim game for the given scenarios.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    final_result_string = ""

    def get_grundy(a):
        """
        Calculates the Grundy value for a pile of size 'a' with a Zero-Move available.
        """
        if a == 0:
            return 0
        if a % 2 == 0:  # a is even and non-zero
            return a - 1
        else:  # a is odd
            return a + 1

    for i, piles in enumerate(scenarios, 1):
        print(f"--- Scenario ({i}): a={piles} ---")
        
        grundy_values = [get_grundy(p) for p in piles]
        
        print("Grundy values g(a) for each pile:")
        equation_parts = []
        for p, g in zip(piles, grundy_values):
            if p == 0:
                 print(f"  g({p}) = 0")
            elif p % 2 == 0:
                 print(f"  g({p}) = {p} - 1 = {g}")
            else:
                 print(f"  g({p}) = {p} + 1 = {g}")
            equation_parts.append(str(g))

        nim_sum = functools.reduce(operator.xor, grundy_values)
        
        nim_sum_str = " ^ ".join(equation_parts)
        print(f"\nNim-sum calculation:")
        print(f"  {nim_sum_str} = {nim_sum}")

        if nim_sum != 0:
            winner = 'A'
            print(f"Result: The Nim-sum is non-zero, so the first player (Alice) wins.\n")
        else:
            winner = 'B'
            print(f"Result: The Nim-sum is zero, so the second player (Bob) wins.\n")
        
        final_result_string += winner

    print("----------------------------")
    print(f"Final combined result string: {final_result_string}")

solve_nim_variant()