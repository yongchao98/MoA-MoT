import functools
import operator

def get_g_value(a):
    """
    Calculates the g-value for a single pile of size 'a'
    with a Zero-Move available.
    """
    if not isinstance(a, int) or a < 0:
        raise ValueError("Pile size must be a non-negative integer.")
    if a == 0:
        return 0
    # For a > 0
    if a % 2 == 1:  # a is odd
        return a + 1
    else:  # a is even
        return a - 1

def solve_and_print_scenario(scenario_num, piles):
    """
    Solves a single scenario, prints the detailed calculation,
    and returns the winner's initial ('A' or 'B').
    """
    # Calculate the g-value for each pile
    g_values = [get_g_value(p) for p in piles]
    
    # Calculate the Nim-sum (XOR sum of all g-values)
    nim_sum = functools.reduce(operator.xor, g_values)

    # Determine the winner
    if nim_sum != 0:
        winner = "Alice"
        winner_char = "A"
    else:
        winner = "Bob"
        winner_char = "B"
        
    # Format the output strings for printing
    pile_str = str(piles)
    g_calculation_str = " XOR ".join(map(str, g_values))
    equation = f"g({piles[0]})"
    for p in piles[1:]:
        equation += f" XOR g({p})"
    
    # Print the detailed results for the scenario
    print(f"({scenario_num}) For piles a = {pile_str}:")
    print(f"    The g-values are: {g_values}")
    print(f"    The Nim-sum is calculated as: {equation} = {g_calculation_str} = {nim_sum}")
    print(f"    The result is {nim_sum}. Since this is {'not zero' if winner == 'Alice' else 'zero'}, {winner} will win.")
    print("-" * 30)

    return winner_char

# --- Main execution ---
scenarios = [
    [12, 12],
    [8, 15, 7],
    [7, 16, 8],
    [12, 19, 21, 10],
    [16, 25, 58, 22, 60]
]

final_result_string = ""
for i, p in enumerate(scenarios, 1):
    final_result_string += solve_and_print_scenario(i, p)

print(f"\nThe final combined result string is: {final_result_string}")
