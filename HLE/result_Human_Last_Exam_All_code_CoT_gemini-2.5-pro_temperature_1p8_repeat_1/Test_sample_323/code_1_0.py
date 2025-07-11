import functools
import operator

def solve_nim_scenarios():
    """
    Analyzes five scenarios of a modified Nim game and determines the winner for each.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    results = []

    print("Analyzing Nim game with Zero-Moves. 'A' for Alice wins, 'B' for Bob wins.\n")

    for i, scenario in enumerate(scenarios):
        a = scenario['a']
        n = scenario['n']
        
        print(f"--- Scenario {i+1}: n={n}, a={a} ---")
        
        # Calculate and print nim_sum
        nim_sum_str = " XOR ".join(map(str, a))
        nim_sum = functools.reduce(operator.xor, a)
        
        print(f"The equation for the Nim-sum is: {nim_sum_str} = {nim_sum}")
        
        winner = ''
        
        # Case 1: Nim-sum is 0
        if nim_sum == 0:
            print(f"The Nim-sum is 0. The winner is determined by the parity of n ({n}).")
            if n % 2 != 0:
                print("Since n is odd, the first player (Alice) will use the last available Zero-Move, forcing Bob into a losing position.")
                winner = 'A'
            else:
                print("Since n is even, the second player (Bob) will use the last available Zero-Move, forcing Alice into a losing position.")
                winner = 'B'
        
        # Case 2: Nim-sum is not 0
        else:
            print(f"The Nim-sum is {nim_sum} (non-zero).")
            # A winning player must move to a state that is a losing position for the opponent.
            # A losing position is defined as having a new Nim-sum of 0 and an even number of piles.
            
            if n % 2 != 0:
                print(f"n ({n}) is odd. To win, Alice must create a state with Nim-sum 0 and an even number of piles.")
                print(f"This requires removing a pile entirely. A pile must exist with size equal to the Nim-sum ({nim_sum}).")
                if nim_sum in a:
                    print(f"A pile of size {nim_sum} exists. Alice removes it, which is a winning move.")
                    winner = 'A'
                else:
                    print(f"No pile of size {nim_sum} exists. Alice cannot force a win, so Bob wins.")
                    winner = 'B'
            else: # n is even
                print(f"n ({n}) is even. To win, Alice must create a state with Nim-sum 0 and an even number of piles.")
                print("This requires modifying a pile without emptying it.")
                
                can_win_without_emptying = False
                # Check if there is a move for Alice to win
                for pile_size in a:
                    if (pile_size ^ nim_sum) < pile_size: # Condition for a valid standard Nim winning move
                        if (pile_size ^ nim_sum) > 0: # Check if the move does not empty the pile
                            can_win_without_emptying = True
                            print(f"Winning move found: Change the pile of size {pile_size} to {pile_size ^ nim_sum}, which is non-zero.")
                            break
                
                if can_win_without_emptying:
                    print("This creates a losing position for Bob (Nim-sum=0, n'=even). Alice wins.")
                    winner = 'A'
                else:
                    print("All moves that make the Nim-sum 0 would empty a pile, creating a winning position for Bob.")
                    winner = 'B'

        print(f"Result for Scenario {i+1}: {winner}")
        results.append(winner)
        print("-" * (20 + len(str(i+1)) + len(str(n)) + len(str(a))))
        print()

    final_answer = "".join(results)
    print(f"The combined result string for all scenarios is: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_nim_scenarios()