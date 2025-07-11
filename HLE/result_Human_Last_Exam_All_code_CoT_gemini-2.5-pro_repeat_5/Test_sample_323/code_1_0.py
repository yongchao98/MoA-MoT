import functools
import operator

def solve_nim_with_zero_moves():
    """
    Solves the modified Nim game for a list of scenarios and prints the results.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]},
    ]
    
    final_result_string = ""
    
    for i, case in enumerate(scenarios, 1):
        n = case['n']
        a = case['a']
        
        # Calculate the Nim-sum of the piles
        nim_sum = functools.reduce(operator.xor, a)
        
        # Create the equation string as requested
        equation_str = " XOR ".join(map(str, a))
        
        # Determine winning conditions
        nim_sum_is_zero = (nim_sum == 0)
        n_is_even = (n % 2 == 0)
        
        # Alice wins if (nim_sum_is_zero) is not equal to (n_is_even)
        if nim_sum_is_zero != n_is_even:
            winner = "Alice"
            final_result_string += 'A'
        else:
            winner = "Bob"
            final_result_string += 'B'
            
        # Print the detailed analysis for the current case
        print(f"({i}) n={n}, a={a}")
        print(f"Equation: {equation_str} = {nim_sum}")
        
        nim_sum_desc = "0" if nim_sum_is_zero else "not 0"
        n_parity_desc = "even" if n_is_even else "odd"
        print(f"Analysis: Nim-sum is {nim_sum_desc} and n is {n_parity_desc}.")
        print(f"Result: {winner} wins.\n")
        
    print(f"Final Answer String: {final_result_string}")

solve_nim_with_zero_moves()