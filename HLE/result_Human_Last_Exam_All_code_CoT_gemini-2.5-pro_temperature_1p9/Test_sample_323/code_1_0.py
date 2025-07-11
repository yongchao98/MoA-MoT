import functools
import operator

def analyze_nim_game(case_num, a):
    """
    Analyzes a single scenario of the modified Nim game, prints the reasoning,
    and returns the winner ('A' for Alice, 'B' for Bob).
    """
    n = len(a)
    nim_sum = functools.reduce(operator.xor, a)
    
    # Print the initial state and Nim-sum calculation
    print(f"({case_num}) n={n}, a={a}")
    equation_parts = [str(x) for x in a]
    equation_str = " XOR ".join(equation_parts)
    print(f"The Nim-sum equation is: {equation_str} = {nim_sum}")

    winner = ''
    # Case 1: The initial Nim-sum is 0.
    if nim_sum == 0:
        # A standard Nim move would give the opponent a winning position.
        # The game becomes a duel of using "Zero-Moves" to pass the turn.
        # With `n` initial Zero-Moves, the first player (Alice) wins if n is odd.
        if n % 2 == 1:
            winner = 'A'  # Alice wins
        else:
            winner = 'B'  # Bob wins
            
    # Case 2: The initial Nim-sum is not 0.
    else:
        # Alice has a standard winning Nim move. She will play to make the Nim-sum 0.
        # She wants to leave Bob in a state where he will lose the subsequent Zero-Move duel.
        # The first player in a Zero-Move duel on `k` items loses if `k` is even.
        # So Alice's goal is to make the number of piles `k'` after her move even.

        if n % 2 == 0:
            # Current number of piles 'n' is even.
            # Alice wants to leave 'n' piles (an even number) for Bob.
            # To do this, her move (p -> p XOR nim_sum) must not empty the pile.
            # This is possible if she can choose a pile `p` such that `p != nim_sum`.
            d = nim_sum.bit_length() - 1
            choosable_piles = [p for p in a if (p >> d) & 1]
            
            # If Alice can find a choosable pile that is not equal to the nim_sum, she wins.
            if any(p != nim_sum for p in choosable_piles):
                winner = 'A'
            else:
                # All choosable piles are equal to the nim_sum. Her only move empties a pile,
                # leaving n-1 piles (an odd number). Bob wins the ZM duel.
                winner = 'B'
        else: # n is odd
            # Current number of piles 'n' is odd.
            # Alice wants to leave 'n-1' piles (an even number) for Bob.
            # To do this, her move must empty the pile, which requires choosing a pile `p` equal to the nim_sum.
            if nim_sum in a:
                winner = 'A'
            else:
                # No such pile exists. Any move she makes leaves 'n' piles (an odd number). Bob wins the ZM duel.
                winner = 'B'
    
    return winner

def main():
    """
    Main function to run all scenarios and print the final result.
    """
    scenarios = [
        {'n': 2, 'a': [12, 12]},
        {'n': 3, 'a': [8, 15, 7]},
        {'n': 3, 'a': [7, 16, 8]},
        {'n': 4, 'a': [12, 19, 21, 10]},
        {'n': 5, 'a': [16, 25, 58, 22, 60]},
    ]
    
    final_result_string = ""
    for i, data in enumerate(scenarios):
        winner = analyze_nim_game(i + 1, data['a'])
        final_result_string += winner
        print(f"Result for ({i+1}): {'Alice' if winner == 'A' else 'Bob'} will win.")
        print("-" * 30)

    print(f"\nThe final concatenated result string is: {final_result_string}")

if __name__ == '__main__':
    main()