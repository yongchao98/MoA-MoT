def solve_mtg_puzzle():
    """
    This function determines the optimal attack and the results of combat.
    """
    
    # Step 1: Identify the creatures to attack with for the optimal play.
    # Since winning this turn is impossible, the optimal play is to attack with all creatures
    # to clear the opponent's board as much as possible, ensuring survival and a win next turn.
    attacking_creatures = [1, 2, 3, 4]
    
    # Step 2: Identify the creatures that will die if the opponent blocks optimally.
    # The opponent's optimal blocks will be:
    # (8) blocks (4) -> (8) dies.
    # (7) blocks (2) -> (2) dies.
    # (6) blocks (1) -> (6) dies.
    # (5) blocks (3) -> neither dies.
    # So, creatures (2), (6), and (8) will die.
    dying_creatures = [2, 6, 8]
    
    # Step 3: Format the output as requested, with numbers in increasing order.
    attacking_creatures.sort()
    dying_creatures.sort()
    
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
    dying_str = ", ".join(f"({c})" for c in dying_creatures)
    
    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()