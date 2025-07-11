def solve_mtg_puzzle():
    """
    This function analyzes the MTG puzzle and determines the optimal attack and resulting casualties.

    The key assumptions are:
    1. Axegrinder Giant (1) has Trample.
    2. Wind Drake (4) has Flying.
    3. Skywinder Drake (8) has Flying and can only block flyers.
    
    The optimal play is to attack with all creatures. The opponent must block to survive,
    leading to a specific set of trades that allows one of your creatures to deal lethal damage.
    """
    
    # Define creatures that should attack to win this turn.
    # Attacking with all creatures is the only way to force a win.
    attacking_creatures = [1, 2, 3, 4]
    
    # Determine which creatures die based on optimal blocking by the opponent.
    # (4) is blocked by (8) -> (4) dies.
    # (1) is double-blocked by (6) and (7) to prevent trample -> (1), (6), (7) die.
    # (2) is blocked by (5) -> (2) and (5) die.
    # (3) is unblocked for lethal damage.
    dying_creatures = [1, 2, 4, 5, 6, 7]
    
    # Sort the lists in increasing order as required.
    attacking_creatures.sort()
    dying_creatures.sort()
    
    # Format the output string as per the user's request.
    attackers_str = ", ".join([f"({c})" for c in attacking_creatures])
    dying_str = ", ".join([f"({c})" for c in dying_creatures])
    
    final_answer = f"{attackers_str}; {dying_str}"
    
    print("Based on the optimal play to win this turn:")
    print(f"Creatures to attack with: {attackers_str}")
    print(f"Creatures that will die: {dying_str}")
    print("\nFinal Answer String:")
    print(final_answer)

solve_mtg_puzzle()