import functools
import operator

def solve_nim_variant():
    """
    Determines the winner for several scenarios of a modified Nim game.
    """
    scenarios = [
        {'id': 1, 'a': [12, 12]},
        {'id': 2, 'a': [8, 15, 7]},
        {'id': 3, 'a': [7, 16, 8]},
        {'id': 4, 'a': [12, 19, 21, 10]},
        {'id': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    final_answer_string = ""
    
    for scenario in scenarios:
        case_num = scenario['id']
        a = scenario['a']
        n = len(a)
        
        # Calculate Nim-sum (s)
        if not a:
            s = 0
        else:
            s = functools.reduce(operator.xor, a)
            
        print(f"Scenario ({case_num}): a = {a}")
        print(f"Number of piles, n = {n}")
        
        # Output the equation for the Nim-sum
        equation = ' ^ '.join(map(str, a)) + f" = {s}"
        print(f"Nim-sum, S = {equation}")

        winner = ''
        
        if s == 0:
            # If S is 0, the winner depends on the parity of n.
            if n % 2 == 0:
                winner = 'A'
                print("Result: Alice wins. (S=0 and n is even)")
            else:
                winner = 'B'
                print("Result: Bob wins. (S=0 and n is odd)")
        else: # s != 0
            if n % 2 == 1: # n is odd, Alice wants n to remain odd after her move.
                # Alice wins if she can make a move without emptying a pile.
                # This requires finding a movable pile a_i (where (a_i^s) < a_i)
                # that is not equal to s.
                can_make_non_emptying_move = False
                for pile_size in a:
                    if (pile_size ^ s) < pile_size:
                        if pile_size != s:
                            can_make_non_emptying_move = True
                            break
                if can_make_non_emptying_move:
                    winner = 'A'
                    print("Result: Alice wins. (S!=0, n is odd, and a non-emptying move is possible)")
                else:
                    winner = 'B'
                    print("Result: Bob wins. (S!=0, n is odd, but all moves must empty a pile)")
            else: # n is even, Alice wants n to become odd after her move.
                # Alice wins if she can empty a pile.
                # This requires finding a pile equal to s.
                if s in a:
                    winner = 'A'
                    print(f"Result: Alice wins. (S!=0, n is even, and a pile of size {s} exists to be removed)")
                else:
                    winner = 'B'
                    print("Result: Bob wins. (S!=0, n is even, but no pile can be removed to make S=0)")

        final_answer_string += winner
        print("-" * 30)

    print(f"The final result string is: {final_answer_string}")

solve_nim_variant()