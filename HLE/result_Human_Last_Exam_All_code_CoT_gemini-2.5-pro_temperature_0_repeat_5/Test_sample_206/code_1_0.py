def solve_mtg_puzzle():
    """
    This function analyzes the MTG board state to determine the optimal attack.
    """

    # Define player and opponent creatures with their attributes
    # Format: {id: {"name": str, "power": int, "toughness": int, "flying": bool}}
    my_creatures = {
        1: {"name": "Axegrinder Giant", "power": 6, "toughness": 4, "flying": False},
        2: {"name": "Centaur Courser", "power": 3, "toughness": 3, "flying": False},
        3: {"name": "Axebane Beast", "power": 2, "toughness": 5, "flying": False},
        4: {"name": "Wind Drake", "power": 2, "toughness": 2, "flying": True},
    }

    opponent_creatures = {
        5: {"name": "River Bear", "power": 2, "toughness": 2, "flying": False},
        6: {"name": "Grizzly Bears", "power": 2, "toughness": 2, "flying": False},
        7: {"name": "Rusted Sentinel", "power": 3, "toughness": 4, "flying": False},
        8: {"name": "Skywinder Drake", "power": 3, "toughness": 1, "flying": True, "block_flying_only": True},
    }

    # --- Strategic Analysis ---
    # A win this turn is impossible because the opponent has a matching number of blockers
    # for both ground and flying attacks.
    # The optimal play is to attack with all ground creatures to eliminate the opponent's
    # ground defense without losing any of our own creatures. This sets up a guaranteed
    # win on the following turn.

    # 1. Determine the attackers
    attackers = [1, 2, 3]

    # 2. Determine the casualties from the optimal blocks
    # Opponent must block all attackers.
    # (7) Rusted Sentinel (3/4) blocks (1) Axegrinder Giant (6/4) -> Sentinel dies.
    # (5) River Bear (2/2) blocks (2) Centaur Courser (3/3) -> Bear dies.
    # (6) Grizzly Bears (2/2) blocks (3) Axebane Beast (2/5) -> Bears die.
    # All of our attacking creatures survive.
    casualties = [5, 6, 7]

    # 3. Format the output string
    # Sort the lists to ensure increasing order
    attackers.sort()
    casualties.sort()

    # Create the formatted strings for each part of the answer
    attackers_str_parts = []
    for creature_id in attackers:
        attackers_str_parts.append(f"({creature_id})")
    
    casualties_str_parts = []
    for creature_id in casualties:
        casualties_str_parts.append(f"({creature_id})")

    attackers_str = ", ".join(attackers_str_parts)
    casualties_str = ", ".join(casualties_str_parts)

    # Print the final result in the specified format
    print(f"{attackers_str}; {casualties_str}")

solve_mtg_puzzle()