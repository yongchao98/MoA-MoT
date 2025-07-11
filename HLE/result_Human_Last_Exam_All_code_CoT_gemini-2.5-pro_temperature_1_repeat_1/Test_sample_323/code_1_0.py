import functools
import operator

def get_grundy_value(a):
    """
    Calculates the Grundy number (g-value) for a single pile of size 'a'
    with one Zero-Move available.
    """
    if a == 0:
        return 0
    if a % 2 == 1:  # a is odd
        return a + 1
    else:  # a is even and > 0
        return a - 1

def solve_nim_variant():
    """
    Solves the five scenarios of the modified Nim game and determines the winner.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    final_result_str = ""
    
    for i, s in enumerate(scenarios):
        piles = s['a']
        print(f"({i+1}) n={s['n']} a={piles}")
        
        g_values = [get_grundy_value(p) for p in piles]
        
        # Format the equation string
        equation_str = " XOR ".join(map(str, g_values))
        
        # Calculate the Nim-sum
        nim_sum = functools.reduce(operator.xor, g_values)
        
        print(f"The g-values for the piles are: {g_values}")
        print(f"The Nim-sum is calculated as: {equation_str} = {nim_sum}")
        
        if nim_sum == 0:
            winner = 'B'
            print("Result: The Nim-sum is 0, so the second player (Bob) will win.\n")
        else:
            winner = 'A'
            print(f"Result: The Nim-sum is {nim_sum} (non-zero), so the first player (Alice) will win.\n")
            
        final_result_str += winner
        
    print(f"The combined result string is: {final_result_str}")
    print(f"<<<{final_result_str}>>>")

solve_nim_variant()