import functools
import operator

def solve_nim_variant():
    """
    Solves the Nim variant problem for the 5 given scenarios,
    printing the logic and the final result.
    """

    def g(a):
        """
        Calculates the Grundy value for a single pile of size 'a'
        in this specific Nim variant.
        """
        if a == 0:
            return 0
        if a % 2 == 1:  # a is odd
            return a + 1
        else:  # a is even and positive
            return a - 1

    scenarios = {
        1: [12, 12],
        2: [8, 15, 7],
        3: [7, 16, 8],
        4: [12, 19, 21, 10],
        5: [16, 25, 58, 22, 60]
    }

    final_result_string = ""
    for i, piles in scenarios.items():
        print(f"--- Scenario ({i}) ---")
        print(f"n={len(piles)} a={piles}")
        
        g_values = []
        
        print("Calculating Grundy values (g-values) for each pile:")
        for pile_size in piles:
            g_val = g(pile_size)
            g_values.append(g_val)
            if pile_size > 0 and pile_size % 2 == 1:
                print(f"g({pile_size}) = {pile_size} + 1 = {g_val} (since pile size is odd)")
            elif pile_size > 0 and pile_size % 2 == 0:
                print(f"g({pile_size}) = {pile_size} - 1 = {g_val} (since pile size is even)")
            else:
                 print(f"g(0) = 0")
        
        nim_sum = functools.reduce(operator.xor, g_values)
        
        g_values_str = " XOR ".join(map(str, g_values))
        print(f"\nTotal Nim-sum = {g_values_str} = {nim_sum}")

        if nim_sum == 0:
            winner = 'B'
            print("Result: Nim-sum is 0, so the second player (Bob) will win.\n")
        else:
            winner = 'A'
            print(f"Result: Nim-sum is non-zero ({nim_sum}), so the first player (Alice) will win.\n")
        
        final_result_string += winner

    print("--- Final Answer ---")
    print("The combined result string for all scenarios is:")
    print(final_result_string)

solve_nim_variant()