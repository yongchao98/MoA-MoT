import functools
import operator

def solve_nim_scenarios():
    """
    Solves the modified Nim game for a list of predefined scenarios.
    For each scenario, it prints a detailed analysis and determines the winner.
    Finally, it prints a concatenated string of the results.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    final_result = []

    print("Analyzing the game scenarios based on the derived winning condition:\n")
    print("Alice (first player) wins if either:")
    print("  - The initial Nim-sum is 0 AND the number of piles (n) is odd.")
    print("  - The initial Nim-sum is not 0 AND the number of piles (n) is even.")
    print("-" * 40)

    for i, scenario in enumerate(scenarios):
        n = scenario['n']
        a = scenario['a']
        
        # Calculate the Nim-sum
        nim_sum = functools.reduce(operator.xor, a)
        
        # Determine winning conditions
        nim_sum_is_zero = (nim_sum == 0)
        n_is_odd = (n % 2 != 0)
        
        # Alice wins if (nim_sum_is_zero == n_is_odd)
        if nim_sum_is_zero == n_is_odd:
            winner = 'A'
            winner_name = "Alice"
        else:
            winner = 'B'
            winner_name = "Bob"
            
        final_result.append(winner)

        # Print the detailed analysis for the current case
        print(f"Scenario ({i+1}): n={n}, a={a}")
        equation = ' XOR '.join(map(str, a))
        print(f"1. Nim-sum calculation: {equation} = {nim_sum}")
        print(f"2. Number of piles n={n} is {'odd' if n_is_odd else 'even'}.")
        
        if winner == 'A':
            if nim_sum_is_zero:
                print("3. Result: Nim-sum is 0 and n is odd. Alice wins.")
            else:
                print("3. Result: Nim-sum is not 0 and n is even. Alice wins.")
        else:
            if nim_sum_is_zero:
                print("3. Result: Nim-sum is 0 and n is even. Bob wins.")
            else:
                print("3. Result: Nim-sum is not 0 and n is odd. Bob wins.")
        print("-" * 40)

    print("Final combined result string:")
    print("".join(final_result))

# Execute the function
solve_nim_scenarios()