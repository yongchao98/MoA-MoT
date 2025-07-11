import functools
import operator

def solve_nim_variant():
    """
    Solves the Nim variant for the 5 given scenarios, explains the reasoning for each,
    and prints the final result string.
    """
    # Define the five scenarios to be analyzed.
    scenarios = {
        1: [12, 12],
        2: [8, 15, 7],
        3: [7, 16, 8],
        4: [12, 19, 21, 10],
        5: [16, 25, 58, 22, 60]
    }

    final_result_string = ""

    # Iterate through each scenario, calculate the result, and print the explanation.
    for i, piles in scenarios.items():
        n = len(piles)

        # 1. Calculate the Nim-sum of the stone piles.
        pile_nim_sum = functools.reduce(operator.xor, piles)

        # 2. Calculate the Nim-sum of the available Zero-Moves.
        # This is 1 if n is odd, 0 if n is even.
        zero_move_nim_sum = n % 2

        # 3. Calculate the total Grundy value of the game state.
        total_nim_sum = pile_nim_sum ^ zero_move_nim_sum

        # 4. Determine the winner based on the total Nim-sum.
        if total_nim_sum != 0:
            winner_char = 'A'
            winner_name = "Alice"
        else:
            winner_char = 'B'
            winner_name = "Bob"
            
        final_result_string += winner_char

        # --- Output the analysis for the current scenario ---
        print(f"--- Scenario ({i}) ---")
        print(f"Piles: a={piles}, Number of piles: n={n}")
        
        # Explain and show the pile Nim-sum calculation
        pile_nim_sum_str = " XOR ".join(map(str, piles))
        print(f"\n1. Pile Nim-sum: {pile_nim_sum_str} = {pile_nim_sum}")
        
        # Explain and show the zero-move Nim-sum calculation
        print(f"2. Zero-move Nim-sum (since n={n}): {n} % 2 = {zero_move_nim_sum}")
        
        # Explain and show the final combined calculation
        print("\n3. Total Nim-sum = (Pile Nim-sum) XOR (Zero-move Nim-sum)")
        print(f"   Final Equation: {pile_nim_sum} XOR {zero_move_nim_sum} = {total_nim_sum}")
        
        # Announce the winner
        print(f"\nResult: The total Nim-sum is {total_nim_sum}. As this is {'non-zero' if total_nim_sum != 0 else 'zero'}, the first player, {winner_name}, wins.")
        print("-" * 25 + "\n")

    # After analyzing all cases, print the final consolidated answer string.
    print(f"The final result string for all scenarios is: {final_result_string}")

# Execute the function to get the solution.
solve_nim_variant()