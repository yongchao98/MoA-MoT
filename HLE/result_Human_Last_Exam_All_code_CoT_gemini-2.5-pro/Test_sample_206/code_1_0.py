def solve_mtg_puzzle():
    """
    This function analyzes the MTG board state and determines the optimal attack and its outcome.

    My Board:
    - Life: 2
    - Creatures:
      (1) Axegrinder Giant (6/4, must attack)
      (2) Centaur Courser (3/3)
      (3) Axebane Beast (3/5)
      (4) Wind Drake (2/2, flying)

    Opponent's Board:
    - Life: 3
    - Creatures:
      (5) River Bear (2/2)
      (6) Grizzly Bears (2/2)
      (7) Rusted Sentinel (3/4, defender)
      (8) Skywinder Drake (3/1, flying, can only block fliers)
    """

    # Step 1: Determine the attackers.
    # The Axegrinder Giant (1) must attack. To secure a win and eliminate threats,
    # the optimal play is an all-out attack with all non-defender creatures.
    # This forces the opponent, who is at a low life total of 3, to make desperate blocks.
    attacking_creatures = [1, 2, 4]
    attacking_creatures.sort()

    # Step 2: Determine the creatures that will die based on optimal blocks.
    # The opponent must block to survive.
    # - Opponent's (8) Skywinder Drake (3/1) must block your (4) Wind Drake (2/2). Both die.
    # - Opponent's (6) Grizzly Bears (2/2) blocks your (2) Centaur Courser (3/3). The bears die, the courser survives.
    # - Opponent's (5) River Bear (2/2) and (7) Rusted Sentinel (3/4) must combine to block your (1) Axegrinder Giant (6/4).
    #   The combined 5 power from the blockers is enough to kill the Giant. The Giant's 6 power is enough to kill both blockers.
    #   So, your Giant (1), their Bear (5), and their Sentinel (7) all die.
    
    dead_creatures = [
        # My creatures that die
        1,  # Axegrinder Giant
        4,  # Wind Drake
        # Opponent's creatures that die
        5,  # River Bear
        6,  # Grizzly Bears
        7,  # Rusted Sentinel
        8,  # Skywinder Drake
    ]
    dead_creatures.sort()

    # Step 3: Format the output as requested.
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
    dead_str = ", ".join(f"({c})" for c in dead_creatures)

    print(f"{attackers_str}; {dead_str}")

solve_mtg_puzzle()