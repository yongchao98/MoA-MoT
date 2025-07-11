def solve_domino_game(n):
    """
    Determines the winner of the domino game for a board of size n.
    
    A position n is a losing position for the current player if n % 3 == 1.
    Otherwise, it is a winning position.
    
    The first player has a winning strategy if n is a winning position.
    The second player has a winning strategy if n is a losing position.
    """
    
    print(f"Analyzing board of size n = {n}")
    
    remainder = n % 3
    print(f"{n} mod 3 = {remainder}")
    
    if remainder == 1:
        print("This is a P-position (losing position for the first player).")
        print("The second player has a winning strategy.")
        return "Second player"
    else:
        print("This is an N-position (winning position for the first player).")
        print("The first player has a winning strategy.")
        return "First player"

print("--- Analysis for n = 20 ---")
winner_20 = solve_domino_game(20)
print("-" * 20)

print("--- Analysis for n = 24 ---")
winner_24 = solve_domino_game(24)
print("-" * 20)

print("\nConclusion:")
print(f"When n = 20, the {winner_20} has a winning strategy.")
print(f"When n = 24, the {winner_24} has a winning strategy.")
