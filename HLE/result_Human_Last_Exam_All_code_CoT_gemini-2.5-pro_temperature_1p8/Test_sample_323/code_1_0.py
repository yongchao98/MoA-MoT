import functools
import operator

def solve_nim_variant():
    """
    Analyzes five scenarios of a modified Nim game and determines the winner for each.

    The game is standard Nim with an added "Zero-Move" option, where a player
    can pass their turn, with the constraint that each pile can only have a
    Zero-Move used on it once. The optimal strategy is determined by the
    interplay between the Nim-sum and the parity of the number of available
    Zero-Moves.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]}
    ]
    
    final_result_string = ""

    for i, sc in enumerate(scenarios):
        n = sc['n']
        a = sc['a']
        
        nim_sum = functools.reduce(operator.xor, a)
        
        equation_str = " XOR ".join(map(str, a)) + f" = {nim_sum}"
        
        print(f"Case ({i+1}): n={n}, a={a}")
        print(f"S = {equation_str}")
        
        winner = ''
        
        if nim_sum == 0:
            print(f"The Nim-sum is 0. The winner is determined by the parity of n={n}.")
            if n % 2 == 1:
                winner = 'A'
                print("n is odd, so Alice can ensure she takes the last Zero-Move. Alice wins.")
            else:
                winner = 'B'
                print("n is even, so Bob can ensure he takes the last Zero-Move. Bob wins.")
        else: # nim_sum > 0
            print(f"The Nim-sum is {nim_sum} (non-zero). Alice will make a standard Nim move to a state with Nim-sum 0.")
            if n % 2 == 1: # n is odd
                print(f"n is odd, so Alice wants to leave Bob with an even number of Zero-Moves ({n-1}).")
                print(f"This requires an emptying move, i.e., finding a pile of size {nim_sum} to remove.")
                if nim_sum in a:
                    winner = 'A'
                    print(f"A pile of size {nim_sum} exists. Alice empties it. Alice wins.")
                else:
                    winner = 'B'
                    print(f"No pile of size {nim_sum} exists. Alice must make a non-emptying move, which lets Bob win. Bob wins.")
            else: # n is even
                print(f"n is even, so Alice wants to leave Bob with an even number of Zero-Moves ({n}).")
                print("This requires a non-emptying move.")
                
                # A winning move on a pile `x` to `x ^ nim_sum` is non-emptying if x != nim_sum.
                # Such a move must also be valid: (x ^ nim_sum) < x.
                non_emptying_move_exists = any((x ^ nim_sum) < x and x != nim_sum for x in a)

                if non_emptying_move_exists:
                    winner = 'A'
                    print("A valid non-emptying winning move exists. Alice makes this move. Alice wins.")
                else:
                    winner = 'B'
                    print("Only emptying winning moves exist. Alice must make one, which lets Bob win. Bob wins.")

        final_result_string += winner
        print("-" * 25)

    print(f"Final Answer String: {final_result_string}")

solve_nim_variant()