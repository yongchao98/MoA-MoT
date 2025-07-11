import functools
import operator

def grundy(s):
    """
    Calculates the Grundy number for a pile of size s with a Zero-Move.
    g(0) = 0
    g(s) = s + 1 if s is odd and s > 0
    g(s) = s - 1 if s is even and s > 0
    """
    if s == 0:
        return 0
    if s % 2 == 1: # Odd
        return s + 1
    else: # Even and positive
        return s - 1

def solve_nim_variant(scenarios):
    """
    Solves the modified Nim game for a list of scenarios.
    
    Args:
        scenarios: A list of tuples, where each tuple is (n, a) representing
                   the number of piles and a list of pile sizes.
    """
    result_string = ""
    for i, (n, piles) in enumerate(scenarios, 1):
        print(f"Scenario ({i}): n={n} a={piles}")
        
        # Calculate Grundy numbers for each pile
        g_values = [grundy(p) for p in piles]
        
        # Calculate the Nim-sum
        nim_sum = functools.reduce(operator.xor, g_values)
        
        # Build the equation string
        g_strings = [f"g({p})" for p in piles]
        g_value_strings = [str(g) for g in g_values]
        
        equation_part1 = " ^ ".join(g_strings)
        equation_part2 = " ^ ".join(g_value_strings)
        
        print(f"  Calculation: {equation_part1} = {equation_part2} = {nim_sum}")
        
        # Determine the winner
        if nim_sum == 0:
            winner = "Bob"
            result_string += "B"
        else:
            winner = "Alice"
            result_string += "A"
            
        print(f"  Result: Nim-sum is {nim_sum}, so {winner} wins.\n")

    print(f"Final answer string: {result_string}")
    return result_string

# Define the scenarios from the problem
scenarios = [
    (2, [12, 12]),
    (3, [8, 15, 7]),
    (3, [7, 16, 8]),
    (4, [12, 19, 21, 10]),
    (5, [16, 25, 58, 22, 60])
]

# Run the solver
final_answer = solve_nim_variant(scenarios)

# Final answer in the required format
print(f"\n<<< {final_answer} >>>")
