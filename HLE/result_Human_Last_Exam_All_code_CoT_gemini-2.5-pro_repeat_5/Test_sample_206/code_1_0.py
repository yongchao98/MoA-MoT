def solve_mtg_puzzle():
    """
    Analyzes the MTG board state and determines the optimal attack and its outcome.
    """
    # Define creatures with their attributes: (name, power, toughness, abilities)
    my_creatures = {
        1: ("Axegrinder Giant", 6, 5, []),
        2: ("Centaur Courser", 3, 3, []),
        3: ("Axebane Beast", 2, 5, []),
        4: ("Wind Drake", 2, 2, ["Flying"]),
    }

    opponent_creatures = {
        5: ("River Bear", 4, 3, ["Islandwalk"]),
        6: ("Grizzly Bears", 2, 2, []),
        7: ("Rusted Sentinel", 3, 4, []),
        8: ("Skywinder Drake", 3, 1, ["Flying", "Can't block non-flyers"]),
    }

    # Step 1: Determine the optimal attackers.
    # The goal is to survive the opponent's next turn (lethal unblockable River Bear)
    # and secure a win. Winning this turn is impossible as the opponent can block all attackers.
    # The optimal play is to attack with all creatures to force the opponent to trade their
    # entire board, including the critical River Bear.
    attacking_creatures_ids = sorted(my_creatures.keys())

    # Step 2: Determine which creatures will die based on optimal blocks.
    # The opponent must block to survive and will make the most logical trades.
    # (8) blocks (4) -> Both die (2/2 vs 3/1)
    # (5) blocks (2) -> Both die (3/3 vs 4/3)
    # (7) blocks (1) -> (7) dies (6/5 vs 3/4)
    # (6) blocks (3) -> (6) dies (2/5 vs 2/2)
    
    my_dead_creatures = [2, 4]
    opponent_dead_creatures = [5, 6, 7, 8]
    
    all_dead_creatures_ids = sorted(my_dead_creatures + opponent_dead_creatures)

    # Step 3: Format the output string as requested.
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures_ids)
    dead_str = ", ".join(f"({c})" for c in all_dead_creatures_ids)

    final_answer = f"{attackers_str}; {dead_str}"
    
    print("Based on the optimal play to ensure victory:")
    print(f"Creatures to attack with: {attackers_str}")
    print(f"Creatures that will die in combat: {dead_str}")
    print("\nFinal formatted answer:")
    print(final_answer)

solve_mtg_puzzle()
<<< (1), (2), (3), (4); (2), (4), (5), (6), (7), (8) >>>