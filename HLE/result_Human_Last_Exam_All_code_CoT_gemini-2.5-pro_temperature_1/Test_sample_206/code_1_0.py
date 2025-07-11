def solve_mtg_puzzle():
    """
    This function analyzes the MTG board state to determine the optimal attack and the resulting casualties.
    """

    # --- Game State ---
    # My Creatures: {ID: (Name, Power, Toughness, Abilities)}
    my_creatures = {
        1: ("Axegrinder Giant", 6, 4, []),
        2: ("Centaur Courser", 3, 3, []),
        3: ("Axebane Beast", 2, 5, []),
        4: ("Wind Drake", 2, 2, ["Flying"]),
    }

    # Opponent's Creatures: {ID: (Name, Power, Toughness, Abilities)}
    opponent_creatures = {
        5: ("River Bear", 4, 3, []),
        6: ("Grizzly Bears", 2, 2, []),
        7: ("Rusted Sentinel", 3, 4, []),
        8: ("Skywinder Drake", 2, 2, ["Flying", "Can only block fliers"]),
    }
    
    opponent_life = 3

    # --- Analysis ---
    # The goal is to force a win or create an overwhelmingly superior board state.
    # An all-out attack is suboptimal, as the opponent can take 2 damage from the weakest
    # creature (Axebane Beast) and survive while making favorable blocks elsewhere.
    # The best attack presents only significant threats, forcing the opponent to block all of them.
    
    # Optimal attackers are those that force blocks to prevent lethal damage.
    # We leave the Axebane Beast (2/5) back to deny the opponent an easy choice.
    attacking_creatures_ids = [1, 2, 4]

    # --- Opponent's Optimal Defense Simulation ---
    # The opponent must block all three attackers to survive.
    # They will make blocks that best preserve their life and board.
    
    dying_creatures_ids = []

    # Combat 1: Flying trade
    # (4) Wind Drake (2/2) is blocked by (8) Skywinder Drake (2/2). Both die.
    dying_creatures_ids.extend([4, 8])

    # Combat 2: Lethal threat #1
    # (2) Centaur Courser (3/3) is blocked by (5) River Bear (4/3). Both die.
    dying_creatures_ids.extend([2, 5])

    # Combat 3: Lethal threat #2
    # (1) Axegrinder Giant (6/4) is blocked by (7) Rusted Sentinel (3/4). 
    # The Sentinel dies, but the Giant survives.
    dying_creatures_ids.append(7)

    # --- Format and Print Output ---
    # Sort the lists in increasing order as requested.
    attacking_creatures_ids.sort()
    dying_creatures_ids.sort()

    # Create the formatted strings for the final answer.
    attackers_str = ", ".join([f"({c_id})" for c_id in attacking_creatures_ids])
    casualties_str = ", ".join([f"({c_id})" for c_id in dying_creatures_ids])

    print(f"{attackers_str}; {casualties_str}")

solve_mtg_puzzle()