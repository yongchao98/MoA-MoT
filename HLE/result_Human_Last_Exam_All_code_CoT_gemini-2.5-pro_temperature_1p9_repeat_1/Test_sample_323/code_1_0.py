import functools
import operator

def get_grundy_value(k):
    """
    Calculates the Grundy value (nim-value) for a single pile of size k.
    G(k) = mex({G(0), ..., G(k-1)} U {k})
    This results in the pattern:
    - G(0) = 0
    - G(k) = k - 1 for k > 0 and even
    - G(k) = k + 1 for k > 0 and odd
    """
    if k == 0:
        return 0
    elif k % 2 == 0: # k is even and positive
        return k - 1
    else: # k is odd
        return k + 1

def solve_nim_variant(piles):
    """
    Calculates the Grundy values for all piles and determines the winner.
    """
    grundy_values = [get_grundy_value(p) for p in piles]
    nim_sum = functools.reduce(operator.xor, grundy_values)
    return nim_sum, grundy_values

def main():
    """
    Solves the 5 given scenarios and prints the results.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]

    final_result_string = ""
    for i, pile_config in enumerate(scenarios, 1):
        nim_sum, grundy_vals = solve_nim_variant(pile_config)
        
        print(f"Case ({i}) a={pile_config}:")
        for pile_size, g_val in zip(pile_config, grundy_vals):
            print(f"G({pile_size}) = {g_val}")

        equation_str = " ^ ".join(map(str, grundy_vals))
        print(f"Nim-sum = {equation_str} = {nim_sum}")
        
        if nim_sum == 0:
            winner = "Bob"
            final_result_string += "B"
        else:
            winner = "Alice"
            final_result_string += "A"
        
        print(f"Result: The nim-sum is {nim_sum}. {winner} will win.\n")
    
    print("Final answer string:")
    print(final_result_string)

if __name__ == "__main__":
    main()
