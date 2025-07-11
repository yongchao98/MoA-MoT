def solve_mtg_puzzle():
    """
    This function determines the optimal attack and the resulting creature deaths for the MTG puzzle.

    The logic is as follows:
    1.  The player must win this turn, otherwise the opponent's unblockable River Bear will win them the game.
    2.  To win, the player must deal 3 damage. The 3/3 Centaur Courser is the key.
    3.  The player must force a situation where the Centaur Courser is unblocked.
    4.  By attacking with the mandatory Axegrinder Giant (1), the Centaur Courser (2), and the Wind Drake (4),
        the player presents multiple threats.
    5.  The opponent must use their flying blocker (8) on the Wind Drake (4).
    6.  The opponent is forced to use all their remaining ground blockers (5, 6, 7) on the massive
        Axegrinder Giant (1) to eliminate the biggest threat on the board.
    7.  This leaves the Centaur Courser (2) unblocked, which deals the winning 3 damage.
    8.  We then calculate which creatures die in combat.
    """

    attackers = [1, 2, 4]
    dead_creatures = [1, 4, 5, 6, 7]

    # Sort the lists in increasing order as requested
    attackers.sort()
    dead_creatures.sort()

    # Format the output string
    attackers_str = ", ".join(f"({n})" for n in attackers)
    dead_creatures_str = ", ".join(f"({n})" for n in dead_creatures)

    print(f"{attackers_str}; {dead_creatures_str}")

solve_mtg_puzzle()