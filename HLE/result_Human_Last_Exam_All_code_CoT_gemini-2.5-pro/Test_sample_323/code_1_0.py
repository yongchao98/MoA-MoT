import functools
import operator

def solve_nim_variant():
    """
    Solves the five scenarios of the modified Nim game and prints the reasoning for each.
    """
    scenarios = [
        (2, [12, 12]),
        (3, [8, 15, 7]),
        (3, [7, 16, 8]),
        (4, [12, 19, 21, 10]),
        (5, [16, 25, 58, 22, 60]),
    ]

    final_result_string = ""
    case_num = 1

    for n, a in scenarios:
        print(f"--- Scenario ({case_num}) ---")
        print(f"Piles a = {a}")
        
        # Calculate Nim-sum and print the equation
        nim_sum = functools.reduce(operator.xor, a)
        equation = " ^ ".join(map(str, a))
        print(f"Number of piles n = {n}")
        print(f"Nim-sum S = {equation} = {nim_sum}")

        winner = ''
        if nim_sum == 0:
            print("The Nim-sum is 0. The winner is determined by the parity of n.")
            if n % 2 == 0:
                winner = 'B'
                print(f"Since n={n} is even, Bob will use the last Zero-Move, forcing Alice to lose. Bob wins.")
            else:
                winner = 'A'
                print(f"Since n={n} is odd, Alice will use the last Zero-Move, forcing Bob to lose. Alice wins.")
        else:
            print("The Nim-sum is not 0. Alice will make a standard move to a state with Nim-sum 0.")
            if n % 2 != 0:
                print(f"n={n} is odd. Alice wants to leave an even number of piles (n-1) for the Zero-Move phase.")
                print("To do this, she must empty a pile, which is possible only if S is one of the piles.")
                if nim_sum in a:
                    winner = 'A'
                    print(f"S={nim_sum} is in the list of piles. Alice can empty this pile and will win.")
                else:
                    winner = 'B'
                    print(f"S={nim_sum} is not in the list of piles. Alice cannot empty a pile and will lose.")
            else: # n is even
                print(f"n={n} is even. Alice wants to leave an even number of piles (n) for the Zero-Move phase.")
                print("To do this, she must make a move that does not empty a pile.")
                can_move_without_emptying = False
                for x in a:
                    if (x ^ nim_sum) < x:  # Check if a move on pile x is valid
                        if x != nim_sum:   # Check if this move empties the pile
                            can_move_without_emptying = True
                            break
                if can_move_without_emptying:
                    winner = 'A'
                    print("Alice can make a move without emptying a pile. She will win.")
                else:
                    winner = 'B'
                    print("All of Alice's valid moves would empty a pile. She will lose.")

        final_result_string += winner
        print(f"Result for scenario ({case_num}): {winner}\n")
        case_num += 1
        
    print("-------------------------")
    print(f"Final combined result string: {final_result_string}")


solve_nim_variant()