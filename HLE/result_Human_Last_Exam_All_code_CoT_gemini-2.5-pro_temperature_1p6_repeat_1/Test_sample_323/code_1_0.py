import functools
import operator

def solve_nim_variant():
    """
    Solves the modified Nim game for a given set of scenarios.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60]
    ]
    
    final_result_string = ""
    
    for i, a in enumerate(scenarios):
        n = len(a)
        nim_sum = functools.reduce(operator.xor, a)
        
        # Determine winner
        alice_wins = False
        if nim_sum != 0:
            if n % 2 == 0:
                alice_wins = True
        else: # nim_sum == 0
            if n % 2 != 0:
                alice_wins = True
        
        winner = 'A' if alice_wins else 'B'
        final_result_string += winner
        
        # Print detailed explanation for each case
        print(f"({i+1}) n={n} a={a}")
        
        # Print the equation
        equation_str = " XOR ".join(map(str, a))
        print(f"Equation: {equation_str} = {nim_sum}")
        
        # Print reasoning
        n_parity = "even" if n % 2 == 0 else "odd"
        print(f"Nim-sum is {nim_sum} and n={n} is {n_parity}.")
        if alice_wins:
            if nim_sum != 0:
                print(f"Rule: S != 0 and n is even -> Alice wins.")
            else:
                print(f"Rule: S == 0 and n is odd -> Alice wins.")
            print(f"Winner: Alice\n")
        else:
            if nim_sum != 0:
                print(f"Rule: S != 0 and n is odd -> Bob wins.")
            else:
                print(f"Rule: S == 0 and n is even -> Bob wins.")
            print(f"Winner: Bob\n")
            
    print(f"Final result string: {final_result_string}")
    print(f"<<<{final_result_string}>>>")

solve_nim_variant()