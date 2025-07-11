import functools
import operator

def solve_nim_variant():
    """
    Solves the modified Nim game for five specific scenarios and prints the reasoning.
    """
    scenarios = [
        [12, 12],
        [8, 15, 7],
        [7, 16, 8],
        [12, 19, 21, 10],
        [16, 25, 58, 22, 60],
    ]

    final_answer_string = ""

    for i, a in enumerate(scenarios):
        n = len(a)
        nim_sum = functools.reduce(operator.xor, a)

        # Build the equation string for printing
        equation = " XOR ".join(map(str, a)) + f" = {nim_sum}"

        print(f"--- Case ({i+1}): n={n}, a={a} ---")
        print(f"Nim-sum: {equation}")

        winner = ''
        if nim_sum == 0:
            print("The Nim-sum is 0. The game becomes a pure 'passing' game.")
            if n % 2 == 0:
                print(f"The number of piles n={n} is even.")
                print("The second player (Bob) can take the last of the n zero-moves, forcing Alice into a losing position.")
                winner = 'B'
            else:  # n is odd
                print(f"The number of piles n={n} is odd.")
                print("The first player (Alice) can take the last of the n zero-moves, forcing Bob into a losing position.")
                winner = 'A'
        else:  # nim_sum != 0
            print(f"The Nim-sum is non-zero ({nim_sum}). Alice must make a move to a state with Nim-sum 0.")
            print("The winner is then decided by the parity of piles left for the subsequent 'passing game'. Alice wins if this number is even.")

            if n % 2 != 0:  # n is odd
                print(f"The number of piles n={n} is odd. Alice needs to leave an even number of piles (n-1).")
                print("This requires emptying a pile of size equal to the Nim-sum.")
                if nim_sum in a:
                    print(f"A pile of size {nim_sum} exists in {a}. Alice removes it.")
                    winner = 'A'
                else:
                    print(f"No pile of size {nim_sum} exists. Alice must leave n={n} (odd) piles, so Bob will win the passing game.")
                    winner = 'B'
            else:  # n is even
                print(f"The number of piles n={n} is even. Alice needs to leave an even number of piles (n).")
                print("This requires making a move that does not empty the pile.")
                
                can_win_without_emptying = False
                # To win, Alice needs to find a pile `p` such that `p XOR nim_sum < p` (a valid Nim move)
                # AND `p != nim_sum` (so the pile is not emptied).
                for pile_size in a:
                    if (pile_size ^ nim_sum) < pile_size and pile_size != nim_sum:
                        can_win_without_emptying = True
                        break
                
                if can_win_without_emptying:
                    print("Alice can make a move that leaves n (even) piles, so she wins.")
                    winner = 'A'
                else:
                    print("All of Alice's standard winning moves would empty a pile, leaving n-1 (odd) piles, so Bob will win the passing game.")
                    winner = 'B'

        print(f"Conclusion: {'Alice' if winner == 'A' else 'Bob'} will win.")
        final_answer_string += winner

    print("\n====================")
    print(f"Final Answer String: {final_answer_string}")
    print("====================")


solve_nim_variant()