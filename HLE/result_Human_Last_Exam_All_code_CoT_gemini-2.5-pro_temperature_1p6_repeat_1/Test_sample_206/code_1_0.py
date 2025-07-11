def solve_mtg_puzzle():
    """
    Solves the MTG combat puzzle.

    This function determines the optimal attack strategy and identifies which creatures
    will die during combat, assuming optimal play from the opponent.
    
    The analysis concludes that winning is impossible with the opponent at 3 life.
    The "optimal" attack is therefore the one that deals the most possible damage (2)
    and removes the maximum number of enemy blockers. This is achieved by attacking
    with all three ground creatures, overwhelming the opponent's two ground blockers.
    """

    # Player's creatures
    my_creatures = {
        1: {"name": "Axegrinder Giant", "power": 6, "toughness": 4},
        2: {"name": "Centaur Courser", "power": 3, "toughness": 3},
        3: {"name": "Axebane Beast", "power": 2, "toughness": 5},
        4: {"name": "Wind Drake", "power": 2, "toughness": 2, "flying": True}
    }

    # Opponent's creatures
    opponent_creatures = {
        5: {"name": "River Bear", "power": 4, "toughness": 3, "can_block": False},
        6: {"name": "Grizzly Bears", "power": 2, "toughness": 2},
        7: {"name": "Rusted Sentinel", "power": 3, "toughness": 4},
        8: {"name": "Skywinder Drake", "power": 2, "toughness": 1, "flying": True}
    }

    # The optimal attack to deal maximum damage is with the three ground creatures.
    attacking_creatures_nums = [1, 2, 3]

    # Opponent optimally blocks the two strongest attackers.
    # Blocker (7) Rusted Sentinel blocks attacker (1) Axegrinder Giant.
    # Blocker (6) Grizzly Bears blocks attacker (2) Centaur Courser.
    # Attacker (3) Axebane Beast is unblocked.
    
    dead_creatures_nums = []

    # Combat: (1) vs (7)
    giant_survives = my_creatures[1]["toughness"] > opponent_creatures[7]["power"]
    sentinel_dies = opponent_creatures[7]["toughness"] <= my_creatures[1]["power"]
    if sentinel_dies:
        dead_creatures_nums.append(7)

    # Combat: (2) vs (6)
    courser_survives = my_creatures[2]["toughness"] > opponent_creatures[6]["power"]
    bears_die = opponent_creatures[6]["toughness"] <= my_creatures[2]["power"]
    if bears_die:
        dead_creatures_nums.append(6)

    # Sort the lists in increasing order
    attacking_creatures_nums.sort()
    dead_creatures_nums.sort()

    # Format the output string
    attackers_str = ", ".join(map(str, attacking_creatures_nums))
    dead_str = ", ".join(map(str, dead_creatures_nums))
    
    # In the final output, each number must be shown
    final_answer = f"({attackers_str}); ({dead_str})"
    print(final_answer)

solve_mtg_puzzle()