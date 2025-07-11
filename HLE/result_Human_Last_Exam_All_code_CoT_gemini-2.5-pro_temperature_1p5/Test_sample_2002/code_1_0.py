import sys

# It's possible to reach high recursion depths, so increase the limit.
sys.setrecursionlimit(2000)

memo = {}

def is_winning(n, m):
    """
    Determines if the current position (n, m) is a winning one.
    A position is winning if there's a move to a losing position for the opponent.
    
    Args:
        n (int): The current number of tokens in the pile.
        m (int): The number of tokens taken by the previous player.
        
    Returns:
        bool: True if the position is winning, False if it's losing.
    """
    if (n, m) in memo:
        return memo[(n, m)]
    
    # Base case: If there are no tokens left, the current player cannot move and thus loses.
    if n == 0:
        return False
        
    # The current player can take x tokens. The maximum they can take is
    # limited by the number of remaining tokens (n) and twice the previous move (2*m).
    limit = min(n, 2 * m)
    
    # Iterate through all possible moves.
    for x in range(1, limit + 1):
        # The next state for the opponent will be (n - x, x).
        # If this position is a losing one for the opponent, then the current
        # position is a winning one.
        if not is_winning(n - x, x):
            memo[(n, m)] = True
            return True
            
    # If all possible moves lead to winning positions for the opponent,
    # then the current position is a losing one.
    memo[(n, m)] = False
    return False

def find_p2_winning_starts(max_T):
    """
    Finds the initial numbers of tokens T for which the second player has a
    winning strategy.
    
    Args:
        max_T (int): The maximum value of T to check.
        
    Returns:
        list: A list of T values for which the second player wins.
    """
    p2_win_Ts = []

    for T in range(1, max_T + 1):
        # A value T is a P2-win if for all of P1's first moves,
        # P2 is left in a winning position.
        
        # P1's first move, x1, can be any number from 1 to T-1.
        # If T=1, P1 has no valid moves (1 <= x1 < 1 is impossible), so P1 loses by default.
        if T == 1:
            p2_win_Ts.append(1)
            continue
        
        is_p1_win_found = False
        for x1 in range(1, T):
            # The state for P2 is (T - x1) tokens, and P1's move was x1.
            # If this position is NOT winning for P2, then P1 has found a way to win.
            if not is_winning(T - x1, x1):
                is_p1_win_found = True
                break
        
        # If after checking all of P1's moves, no path to a P1 win was found,
        # then T is a P2-win number.
        if not is_p1_win_found:
            p2_win_Ts.append(T)
            
    return p2_win_Ts

def main():
    """
    Main function to run the analysis and print the conclusion.
    """
    max_tokens_to_check = 60
    print(f"Analyzing the game for T up to {max_tokens_to_check}...")
    
    winning_Ts = find_p2_winning_starts(max_tokens_to_check)
    
    print("\nThe second player has a winning strategy if the initial number of tokens T is one of:")
    print(winning_Ts)
    
    # For confirmation, we generate Fibonacci numbers and compare.
    fibs = [1, 2]
    while fibs[-1] <= max_tokens_to_check:
        fibs.append(fibs[-1] + fibs[-2])
    if fibs[-1] > max_tokens_to_check:
        fibs.pop()
    
    print("\nThis sequence of numbers corresponds to the Fibonacci numbers:")
    print(fibs)
    
    if winning_Ts == fibs:
        print("\nConclusion: The second player has a winning strategy if, and only if,")
        print("the initial number of tokens T is a Fibonacci number.")
    else:
        print("\nConclusion: The sequence does not perfectly match Fibonacci numbers.")

if __name__ == "__main__":
    main()
