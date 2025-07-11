def solve_mtg_puzzle():
    """
    Analyzes the MTG board state to determine the optimal attack and the outcome.

    The primary tactical consideration is that the opponent's River Bear has
    Islandwalk, making it unblockable since we control Islands. With our life total
    at 2, this creature represents a lethal threat on the opponent's next turn.

    Therefore, the optimal play is an "alpha strike" (attacking with all creatures)
    to force the opponent to block with their entire board, including the
    River Bear, to prevent losing to combat damage this turn.
    """

    # List of creatures to attack with, determined by the strategy.
    # Attacking with all is necessary to force blocks and remove the threatening River Bear.
    attackers = [1, 2, 3, 4]

    # List of creatures that will die during combat if played optimally.
    # (4) Wind Drake (2/2) is blocked by (8) Skywinder Drake (2/2) -> Both die.
    # (1) Axegrinder Giant (6/4) is blocked by (7) Rusted Sentinel (3/4) -> (7) dies.
    # (2) Centaur Courser (3/3) is blocked by (5) River Bear (2/2) -> (5) dies.
    # (3) Axebane Beast (2/5) is blocked by (6) Grizzly Bears (2/2) -> (6) dies.
    dead_creatures = [4, 5, 6, 7, 8]

    # Sort the lists in increasing order as per the format requirement.
    attackers.sort()
    dead_creatures.sort()

    # Format the output string.
    attackers_str = ", ".join(f"({num})" for num in attackers)
    dead_creatures_str = ", ".join(f"({num})" for num in dead_creatures)

    # Print the final answer in the specified format.
    print(f"{attackers_str}; {dead_creatures_str}")

solve_mtg_puzzle()