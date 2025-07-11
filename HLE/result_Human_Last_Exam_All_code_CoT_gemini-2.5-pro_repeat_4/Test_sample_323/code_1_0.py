import functools
import operator

def g(a):
    """
    Calculates the Grundy value for a single pile of size 'a'
    in the modified Nim game.
    """
    if a == 0:
        return 0
    if a % 2 == 1:  # a is odd
        return a + 1
    else:  # a is even and positive
        return a - 1

def solve_and_print(case_num, piles):
    """
    Solves a single case of the game, prints the step-by-step calculation,
    and returns the winner's initial ('A' for Alice, 'B' for Bob).
    """
    print(f"--- Case ({case_num}): n={len(piles)} a={piles} ---")
    
    # Calculate Grundy values for all piles
    grundy_values = [g(p) for p in piles]
    
    # Print the g-value calculation for each pile
    for i in range(len(piles)):
        print(f"g({piles[i]}) = {grundy_values[i]}")
        
    # Calculate the Nim-sum (XOR sum) of all Grundy values
    nim_sum = functools.reduce(operator.xor, grundy_values)
    
    # Format the equation string to be printed
    equation_str = " XOR ".join(map(str, grundy_values))
    print(f"Nim-sum = {equation_str} = {nim_sum}")
    
    # Determine the winner based on the Nim-sum
    if nim_sum == 0:
        winner = "Bob"
        winner_char = "B"
        print(f"Result: The Nim-sum is 0, so the second player (Bob) will win.")
    else:
        winner = "Alice"
        winner_char = "A"
        print(f"Result: The Nim-sum is non-zero, so the first player (Alice) will win.")
    
    return winner_char

def main():
    """
    Main function to run all test cases and produce the final result string.
    """
    test_cases = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]
    
    final_result_string = ""
    for i, piles in enumerate(test_cases, 1):
        final_result_string += solve_and_print(i, piles)
    
    print("\n" + "="*40)
    print(f"Final combined result string: {final_result_string}")
    print("="*40)


if __name__ == "__main__":
    main()