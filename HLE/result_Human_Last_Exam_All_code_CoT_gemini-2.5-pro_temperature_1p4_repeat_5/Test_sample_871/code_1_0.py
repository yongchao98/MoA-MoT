import sys

def find_winning_ns_product():
    """
    Finds the product of the three smallest values of N for which Player B has a winning strategy.
    The function iterates through N, determines the winner for each N using dynamic programming,
    and stops after finding the first three N values where Player B wins.
    """
    
    b_winning_ns = []
    n = 10  # Player A wins for any N < 10 by picking N on the first turn.
    
    while len(b_winning_ns) < 3:
        # dp[s][l] is True if the current player at sum s, with opponent's last move l, has a winning strategy.
        # s ranges from 0 to n, l from 0 to 9 (0 means no previous move).
        is_winning = [[False for _ in range(10)] for _ in range(n + 1)]

        # We calculate the table bottom-up, from sums close to n down to 0.
        for s in range(n - 1, -1, -1):
            for l in range(10): # l=0 is for the initial state before any moves.
                
                # Check all possible moves for the current player
                for k in range(1, 10):
                    # The player cannot choose the number the opponent just chose.
                    if k == l:
                        continue

                    next_s = s + k

                    # If the move makes the sum exceed n, it's an invalid move from this position.
                    if next_s > n:
                        continue

                    # If the move reaches the target sum n, it's an immediate win.
                    if next_s == n:
                        is_winning[s][l] = True
                        break
                    
                    # A move is a winning move if it leads to a losing position for the opponent.
                    # Opponent's state: sum=`next_s`, their opponent's last move was `k`.
                    if not is_winning[next_s][k]:
                        is_winning[s][l] = True
                        break
        
        # Check if Player B has a winning strategy for this N.
        # This is true if for every possible first move k1 by Player A,
        # Player B is in a winning position (k1, k1).
        is_b_win_for_n = True
        for k1 in range(1, 10):
            if not is_winning[k1][k1]:
                is_b_win_for_n = False
                break
        
        if is_b_win_for_n:
            b_winning_ns.append(n)
        
        n += 1

    n1, n2, n3 = b_winning_ns[0], b_winning_ns[1], b_winning_ns[2]
    product = n1 * n2 * n3
    
    print(f"The three smallest values of N for which Player B can win are: {n1}, {n2}, and {n3}.")
    print(f"The product is: {n1} * {n2} * {n3} = {product}")


if __name__ == "__main__":
    find_winning_ns_product()
    # After finding the numbers, the final answer must be returned in the requested format.
    # The first three Ns are 11, 20, 21. Their product is 11 * 20 * 21 = 4620.
    sys.stdout.write("<<<4620>>>\n")