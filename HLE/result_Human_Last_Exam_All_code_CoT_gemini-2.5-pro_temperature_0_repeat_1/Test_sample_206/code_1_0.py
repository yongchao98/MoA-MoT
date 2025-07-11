def solve_mtg_puzzle():
    """
    This function determines the optimal attack and the outcome of combat.
    """
    
    # List of creatures to attack with, in increasing order of their number.
    attacking_creatures = [1, 2, 4]
    
    # List of all creatures that will die in combat if the opponent blocks optimally.
    # My creatures: (1) Giant, (2) Courser, (4) Drake
    # Opponent's creatures: (5) Bear, (6) Bears, (7) Sentinel, (8) Drake
    dying_creatures = [1, 2, 4, 5, 6, 7, 8]
    
    # Format the output string as per the example.
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
    dying_str = ", ".join(f"({c})" for c in dying_creatures)
    
    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()