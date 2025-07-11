import sys

# It's possible to hit the recursion limit for larger numbers, so we can increase it.
# sys.setrecursionlimit(2000)

memo = {}

def is_losing_position(tokens, prev_move):
    """
    Determines if the current game state is a losing one for the current player.
    A state is defined by (tokens, prev_move).
    A state is losing if every possible move leads to a winning state for the opponent.
    This function uses memoization to store results and avoid re-computation.
    """
    if (tokens, prev_move) in memo:
        return memo[(tokens, prev_move)]

    # If the current player can take all remaining tokens, it's a winning position.
    # The number of tokens they can take is up to 2 * prev_move.
    if tokens <= 2 * prev_move:
        memo[(tokens, prev_move)] = False
        return False

    # A state is a winning position if there is at least one move that leads
    # to a losing position for the opponent.
    # We iterate through all possible moves for the current player.
    max_move = 2 * prev_move
    for move in range(1, max_move + 1):
        if move > tokens:
            break
        
        # The opponent's next state will be (tokens - move, move).
        # If we can make a move that puts the opponent in a losing position...
        if is_losing_position(tokens - move, move):
            # ...then this current state is a winning one for us.
            memo[(tokens, prev_move)] = False
            return False

    # If all of our possible moves lead to winning positions for the opponent,
    # then our current state is a losing one.
    memo[(tokens, prev_move)] = True
    return True

def does_p2_win(T):
    """
    Checks if the second player has a winning strategy for an initial pile of T tokens.
    P2 wins if for every possible first move by P1, P2 is left in a winning state.
    """
    # P1 can take x1 tokens, where 1 <= x1 < T
    for p1_move in range(1, T):
        tokens_left = T - p1_move
        
        # It's now P2's turn. The state is (tokens_left, p1_move).
        # If this state is a losing one for P2, then P1 has found a winning line of play.
        # A winning state for P2 is one that is NOT a losing position.
        if is_losing_position(tokens_left, p1_move):
            # T is a win for P1, not P2.
            return False
    
    # If all of P1's opening moves lead to winning positions for P2, then P2 wins.
    return True

def main():
    """
    Finds and prints all values of T for which the second player has a winning strategy,
    up to a defined limit.
    """
    print("The values of T for which the second player has a winning strategy are:")
    
    # We check for T starting from 2, up to a reasonable limit to demonstrate the pattern.
    limit = 60
    p2_winning_t = []
    
    for T in range(2, limit + 1):
        # Clear memoization cache for each new T to ensure clean calculation.
        memo.clear()
        if does_p2_win(T):
            p2_winning_t.append(T)
            
    # The question asks to output each number in the final equation/sequence.
    for number in p2_winning_t:
        print(number)

    print("\nThese numbers correspond to the Fibonacci sequence.")

if __name__ == "__main__":
    main()
