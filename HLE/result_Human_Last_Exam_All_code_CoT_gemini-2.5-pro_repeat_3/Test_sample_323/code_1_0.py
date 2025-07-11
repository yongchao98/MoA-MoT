import functools
import operator

def get_grundy_value(x):
    """
    Calculates the Grundy value (g-value) for a single pile of size x
    in the modified Nim game.
    """
    if x == 0:
        return 0
    if x % 2 == 1:  # x is odd
        return x + 1
    else:  # x is even and positive
        return x - 1

def solve_nim_variant(piles, case_num):
    """
    Solves a single scenario of the modified Nim game.
    """
    if not piles:
        # No piles, the first player has no move and loses.
        return 'B'

    # 1. Calculate g-value for each pile
    grundy_values = [get_grundy_value(p) for p in piles]

    # 2. Compute the Nim-sum of the g-values
    nim_sum = functools.reduce(operator.xor, grundy_values)

    # 3. Determine the winner
    winner = "Alice" if nim_sum != 0 else "Bob"

    # 4. Format the output string for explanation
    g_piles_str = " XOR ".join([f"g({p})" for p in piles])
    g_values_str = " XOR ".join(map(str, grundy_values))
    
    print(f"({case_num}) a={piles}")
    print(f"    Calculation: {g_piles_str} = {g_values_str} = {nim_sum}")
    print(f"    Result: Nim-sum is {'non-zero' if nim_sum != 0 else 'zero'}, so {winner} will win.")
    
    return winner[0]

def main():
    """
    Main function to solve all provided scenarios and print the final result string.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]
    
    final_result_string = ""
    for i, p in enumerate(scenarios, 1):
        final_result_string += solve_nim_variant(p, i)
        print("-" * 20)

    print(f"\nFinal answer string: {final_result_string}")


if __name__ == "__main__":
    main()