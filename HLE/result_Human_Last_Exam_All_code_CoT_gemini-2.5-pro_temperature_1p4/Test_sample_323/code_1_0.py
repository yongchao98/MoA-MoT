import functools
import operator

def get_grundy_value(k):
    """
    Calculates the Grundy value for a single pile of size k in the modified Nim game.
    g(0) = 0
    g(k) = k + 1 for k odd
    g(k) = k - 1 for k > 0 and even
    """
    if k == 0:
        return 0
    elif k % 2 == 1:  # k is odd
        return k + 1
    else:  # k is even and positive
        return k - 1

def solve_nim_scenarios():
    """
    Solves the five given Nim game scenarios and prints the results.
    """
    scenarios = {
        1: [12, 12],
        2: [8, 15, 7],
        3: [7, 16, 8],
        4: [12, 19, 21, 10],
        5: [16, 25, 58, 22, 60]
    }

    final_answer = ""

    for i in range(1, 6):
        piles = scenarios[i]
        
        # Calculate Grundy values for each pile
        grundy_values = [get_grundy_value(p) for p in piles]
        
        # Calculate the Nim-sum (XOR sum of all Grundy values)
        nim_sum = functools.reduce(operator.xor, grundy_values)
        
        # Determine the winner
        if nim_sum == 0:
            winner = "Bob"
            final_answer += "B"
        else:
            winner = "Alice"
            final_answer += "A"
            
        # Format the output strings for the equation
        g_piles_str = " XOR ".join([f"g({p})" for p in piles])
        g_vals_str = " XOR ".join(map(str, grundy_values))

        print(f"Scenario ({i}) with piles {piles}:")
        print(f"The winning condition is determined by the Nim-sum of the Grundy values.")
        print(f"Calculation: {g_piles_str}")
        print(f"           = {g_vals_str}")
        print(f"           = {nim_sum}")
        print(f"Result: The Nim-sum is {nim_sum}. {winner} will win.\n")

    print(f"The final result string is: {final_answer}")

solve_nim_scenarios()