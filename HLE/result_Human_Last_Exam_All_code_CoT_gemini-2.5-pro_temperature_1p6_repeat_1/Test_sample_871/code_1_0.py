import sys

# Set a higher recursion limit for safety, as the solution involves deep recursion.
sys.setrecursionlimit(2000)

memo = {}
target_N = 0

def can_win(current_sum, opponent_last_move):
    """
    Determines if the current player can win from the given state.
    A state is defined by the current running total and the opponent's last move.
    Uses memoization to store and retrieve results for previously seen states.

    Returns:
        bool: True if the current position is a winning one, False otherwise.
    """
    global memo, target_N
    
    state = (current_sum, opponent_last_move)
    if state in memo:
        return memo[state]

    # First, check if the current player has any valid move at all.
    # A move is invalid if it makes the sum exceed the target N.
    has_at_least_one_valid_move = False
    for my_move in range(1, 10):
        if my_move == opponent_last_move:
            continue
        if current_sum + my_move <= target_N:
            has_at_least_one_valid_move = True
            break
            
    if not has_at_least_one_valid_move:
        # If no move is possible without exceeding N, this is a losing position.
        memo[state] = False
        return False

    # A player wins if they can make a move that leads to a losing position for the opponent.
    # A position is losing for the opponent if can_win(new_sum, my_move) is False.
    for my_move in range(1, 10):
        if my_move == opponent_last_move:
            continue
        
        new_sum = current_sum + my_move
        
        # We only consider moves that do not exceed the target sum.
        if new_sum <= target_N:
            # If the opponent CANNOT win from the resulting state, it means we've found a winning move.
            if not can_win(new_sum, my_move):
                memo[state] = True
                return True

    # If all valid moves lead to a winning position for the opponent,
    # then the current position is a losing one.
    memo[state] = False
    return False

def is_B_win_scenario(N):
    """
    Checks if a given target N results in a win for Player B, assuming optimal play.
    This is true if Player A (the first player) cannot force a win, regardless of their first move.
    """
    global memo, target_N
    # Reset memoization cache and set the global target_N for the new scenario.
    memo = {}
    target_N = N
    
    # Player B wins if for ALL of Player A's possible first moves 'k', 
    # the resulting state is a winning one for Player B.
    for k in range(1, 10):
        # After A moves 'k', the state is (sum=k, last_move=k). It's B's turn.
        # If B cannot win from this state, it means A found a winning opening move.
        if not can_win(k, k):
            return False
            
    # If Player B can win from all possible starting states created by Player A, then B has a winning strategy.
    return True

def solve():
    """
    Finds the three smallest values of N for which Player B has a winning strategy
    and prints their product.
    """
    b_win_Ns = []
    n = 1
    # We loop until we find the three required values of N.
    while len(b_win_Ns) < 3:
        # Player A wins for N=1 to 9 immediately by picking N. We can skip these checks.
        if n > 9:
            if is_B_win_scenario(n):
                b_win_Ns.append(n)
        n += 1
        
    n1, n2, n3 = b_win_Ns[0], b_win_Ns[1], b_win_Ns[2]
    product = n1 * n2 * n3
    
    print("The three smallest values of N for which Player B wins are:")
    print(f"N1 = {n1}")
    print(f"N2 = {n2}")
    print(f"N3 = {n3}")
    print("\nThe product of these three values is the final answer.")
    print(f"Final Equation: {n1} * {n2} * {n3} = {product}")
    
# Execute the solver function
solve()
<<<8228>>>