import functools
import operator

def solve_nim_scenarios():
    """
    Solves the five modified Nim game scenarios and prints the results.
    """
    scenarios = [
        (2, [12, 12]),
        (3, [8, 15, 7]),
        (3, [7, 16, 8]),
        (4, [12, 19, 21, 10]),
        (5, [16, 25, 58, 22, 60]),
    ]

    result_string = ""
    print("Calculating the winner for each scenario:")
    
    for i, (n, a) in enumerate(scenarios):
        # Calculate the Nim-sum of the piles
        nim_sum = functools.reduce(operator.xor, a)
        
        # Calculate the total Grundy value
        g_total = nim_sum ^ n
        
        # Determine the winner
        if g_total != 0:
            winner = "Alice"
            result_string += "A"
        else:
            winner = "Bob"
            result_string += "B"
            
        # Format the equation string
        equation_str = " XOR ".join(map(str, a))
        
        print(f"\nScenario ({i+1}): n={n}, a={a}")
        print(f"Equation: ({equation_str}) XOR {n}")
        print(f"Result: {nim_sum} XOR {n} = {g_total}. Since the result is not 0, {winner} wins.")

    print("\n" + "="*20)
    print("Final answer string:")
    print(result_string)
    
solve_nim_scenarios()