import sys

# Increase recursion limit for deep calculations, although memoization helps a lot.
sys.setrecursionlimit(2000)

memo = {}

def get_grundy(n):
    """
    Calculates the Grundy number (g-value) for a board of size n using recursion and memoization.
    """
    if n in memo:
        return memo[n]
    if n <= 0:
        return 0
    
    reachable_grundy_values = set()
    
    # Iterate through all possible moves to find the g-values of reachable states.
    
    # Case 1: Placing a 1-tile domino.
    # A 1-tile placed at position i (1-indexed) splits the board into two sub-boards
    # of sizes (i-1) and (n-i).
    for i in range(1, n + 1):
        nim_sum = get_grundy(i - 1) ^ get_grundy(n - i)
        reachable_grundy_values.add(nim_sum)
        
    # Case 2: Placing a 2-tile domino.
    # A 2-tile placed at position i (covering i and i+1) splits the board into
    # sub-boards of sizes (i-1) and (n-(i+1)).
    if n >= 2:
        for i in range(1, n):
            nim_sum = get_grundy(i - 1) ^ get_grundy(n - (i + 1))
            reachable_grundy_values.add(nim_sum)
            
    # The g-value is the "minimum excluded value" (mex) of the set of reachable g-values.
    mex = 0
    while mex in reachable_grundy_values:
        mex += 1
    
    memo[n] = mex
    return mex

def solve():
    """
    Solves the problem for n=20 and n=24.
    """
    n1 = 20
    n2 = 24
    
    # Calculate Grundy numbers up to the maximum n required.
    max_n = max(n1, n2)
    for i in range(1, max_n + 1):
        get_grundy(i)
        
    g_n1 = memo[n1]
    g_n2 = memo[n2]
    
    # Determine the winner for n=20
    print(f"For n = {n1}, the Grundy value g({n1}) is {g_n1}.")
    if g_n1 > 0:
        winner_n1 = "first player"
    else:
        winner_n1 = "second player"
    print(f"A Grundy value greater than 0 means it is a winning position for the first player.")
    print(f"Therefore, when n = {n1}, the {winner_n1} has a winning strategy.")

    print("-" * 20)

    # Determine the winner for n=24
    print(f"For n = {n2}, the Grundy value g({n2}) is {g_n2}.")
    if g_n2 > 0:
        winner_n2 = "first player"
    else:
        winner_n2 = "second player"
    print(f"A Grundy value greater than 0 means it is a winning position for the first player.")
    print(f"Therefore, when n = {n2}, the {winner_n2} has a winning strategy.")

    print("-" * 20)

    # Final conclusion matching the answer choices
    print("Conclusion:")
    print(f"When n = {n1}, the {winner_n1} has a winning strategy.")
    print(f"When n = {n2}, the {winner_n2} has a winning strategy.")

solve()