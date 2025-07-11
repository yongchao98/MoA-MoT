def solve_mtg_puzzle():
    """
    This function analyzes the MTG board state to determine the optimal attack
    and the resulting dead creatures, assuming optimal play from both sides.
    """

    # Player's Creatures (You)
    # (1) Axegrinder Giant (6/5), (2) Centaur Courser (3/3)
    # (3) Axebane Beast (2/5), (4) Wind Drake (2/2, Flying)
    
    # Opponent's Creatures
    # (5) River Bear (4/3), (6) Grizzly Bears (2/2)
    # (7) Rusted Sentinel (3/4), (8) Skywinder Drake (3/1, Flying)

    # Since no immediate win is possible, the best strategic move is an all-out attack
    # to maximize pressure and improve the board state for future turns.
    creatures_to_attack_with = [1, 2, 3, 4]

    # Based on optimal blocking by the opponent to ensure survival and minimize losses:
    # - (4) Wind Drake is blocked by (8) Skywinder Drake. Both die.
    # - (2) Centaur Courser is blocked by (7) Rusted Sentinel. (2) dies, (7) survives.
    # - (1) Axegrinder Giant is blocked by (6) Grizzly Bears. (6) dies, (1) survives.
    # - (3) Axebane Beast is blocked by (5) River Bear. Both survive this combat.
    creatures_that_will_die = [2, 4, 6, 8]

    # Sort the lists in increasing order as required.
    creatures_to_attack_with.sort()
    creatures_that_will_die.sort()

    # Format the output string as per the example: (1), (2); (1), (6)
    attackers_str = ", ".join(f"({c})" for c in creatures_to_attack_with)
    dead_str = ", ".join(f"({c})" for c in creatures_that_will_die)

    # Print the final formatted answer
    print(f"{attackers_str}; {dead_str}")

solve_mtg_puzzle()