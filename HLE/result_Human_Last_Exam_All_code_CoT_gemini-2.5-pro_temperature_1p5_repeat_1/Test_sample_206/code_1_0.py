def solve_mtg_puzzle():
    """
    This script determines the optimal attack and the resulting creatures that die
    for the given Magic: The Gathering board state.
    """

    # Define the board state
    my_creatures = {
        1: {'name': 'Axegrinder Giant', 'power': 6, 'toughness': 4},
        2: {'name': 'Centaur Courser', 'power': 3, 'toughness': 3},
        3: {'name': 'Axebane Beast', 'power': 2, 'toughness': 4},
        4: {'name': 'Wind Drake', 'power': 2, 'toughness': 2},
    }

    opponent_creatures = {
        5: {'name': 'River Bear', 'power': 2, 'toughness': 2},
        6: {'name': 'Grizzly Bears', 'power': 2, 'toughness': 2},
        7: {'name': 'Rusted Sentinel', 'power': 3, 'toughness': 4},
        8: {'name': 'Skywinder Drake', 'power': 3, 'toughness': 1},
    }

    # The optimal play is to attack with all creatures to maximize pressure and
    # force the opponent into a board state from which they cannot recover.
    attacking_creatures_ids = sorted(my_creatures.keys())

    # Based on the optimal blocking strategy from the opponent to survive,
    # we can determine which creatures will die.
    # (4) Wind Drake is blocked by (8) Skywinder Drake. 3 damage to the 2/2 kills it. 2 damage to the 3/1 kills it. Both die.
    # (1) Giant is blocked by (7) Sentinel. 6 damage kills the 3/4 Sentinel. 3 damage is not enough to kill the 6/4 Giant.
    # (2) Courser is blocked by a 2/2 Bear (e.g., 5). 3 damage kills the 2/2. 2 damage doesn't kill the 3/3.
    # (3) Beast is blocked by the other 2/2 Bear (e.g., 6). 2 damage kills the 2/2. 2 damage doesn't kill the 2/4.
    
    # List of creatures that will die in combat.
    creatures_that_die = [
        4,  # Your Wind Drake
        5,  # Opponent's River Bear
        6,  # Opponent's Grizzly Bears
        7,  # Opponent's Rusted Sentinel
        8,  # Opponent's Skywinder Drake
    ]
    
    # Sort the lists and format for the final output.
    attacking_str = ", ".join(f"({i})" for i in sorted(attacking_creatures_ids))
    died_str = ", ".join(f"({i})" for i in sorted(creatures_that_die))

    final_answer = f"{attacking_str}; {died_str}"
    
    print(final_answer)

solve_mtg_puzzle()