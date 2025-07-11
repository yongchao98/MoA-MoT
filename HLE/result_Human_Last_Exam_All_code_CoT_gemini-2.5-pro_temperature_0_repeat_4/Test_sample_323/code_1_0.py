import functools
import operator

def solve_nim_variant():
    """
    Solves a series of Nim game variants with a "Zero-Move" rule.
    For each scenario, it determines the winner and prints the reasoning.
    Finally, it prints a consolidated string of winners.
    """
    scenarios = [
        (2, [12, 12]),
        (3, [8, 15, 7]),
        (3, [7, 16, 8]),
        (4, [12, 19, 21, 10]),
        (5, [16, 25, 58, 22, 60])
    ]
    
    final_result_string = ""

    for i, (n, a) in enumerate(scenarios):
        case_num = i + 1
        
        # Calculate the Nim-sum (XOR sum) of the piles
        nim_sum = functools.reduce(operator.xor, a)
        
        winner = ''
        
        print(f"Case ({case_num}): n={n}, a={a}")
        
        # Print the equation for the Nim-sum
        equation_str = " XOR ".join(map(str, a))
        print(f"Equation: {equation_str} = {nim_sum}")

        if nim_sum == 0:
            # If Nim-sum is 0, the game is decided by the parity of n (number of Zero-Moves)
            if n % 2 == 1:
                winner = 'A'
                print(f"Nim-sum is 0 and n ({n}) is odd.")
            else:
                winner = 'B'
                print(f"Nim-sum is 0 and n ({n}) is even.")
        else: # nim_sum != 0
            if n % 2 == 1:
                # Alice needs to make n even, by emptying a pile.
                # This is possible iff a pile size equals the nim_sum.
                if nim_sum in a:
                    winner = 'A'
                    print(f"Nim-sum is {nim_sum}, n ({n}) is odd, and a pile of size {nim_sum} exists.")
                else:
                    winner = 'B'
                    print(f"Nim-sum is {nim_sum}, n ({n}) is odd, and no pile of size {nim_sum} exists.")
            else: # n % 2 == 0
                # Alice can win by making a non-emptying move unless a special condition is met.
                # The special condition forces any S'=0 move to be an emptying move.
                can_make_non_emptying_move = False
                for x in a:
                    # A move on pile x to make nim_sum 0 requires changing x to (nim_sum ^ x).
                    # The move is valid if (nim_sum ^ x) < x.
                    # The move is non-emptying if (nim_sum ^ x) > 0, i.e., nim_sum != x.
                    if (nim_sum ^ x) < x:
                        if nim_sum != x:
                            can_make_non_emptying_move = True
                            break
                
                if can_make_non_emptying_move:
                    winner = 'A'
                    print(f"Nim-sum is {nim_sum}, n ({n}) is even, and Alice can make a non-emptying move to a state with Nim-sum 0.")
                else:
                    winner = 'B'
                    print(f"Nim-sum is {nim_sum}, n ({n}) is even, but any move to a state with Nim-sum 0 empties a pile.")

        if winner == 'A':
            print("Result: Alice wins.\n")
        else:
            print("Result: Bob wins.\n")
            
        final_result_string += winner

    print("Final answer string:")
    print(final_result_string)

solve_nim_variant()