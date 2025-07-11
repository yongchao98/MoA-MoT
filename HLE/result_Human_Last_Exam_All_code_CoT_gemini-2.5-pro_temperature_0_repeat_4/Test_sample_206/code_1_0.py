def solve_mtg_puzzle():
    """
    This function determines the optimal attack strategy for the given MTG scenario.

    The logic is as follows:
    1. The player is at 2 life, and the opponent has a 2/2 flyer ((8) Skywinder Drake).
       If this flyer is not dealt with, the opponent will win on their next turn.
    2. To deal with the flyer, the player must attack with their own flyer ((4) Wind Drake),
       forcing a block and trade with the opponent's flyer.
    3. To maximize board advantage and set up a win for the next turn, the player should
       attack with all other creatures as well.
    4. The opponent, to survive the turn, will block all incoming attackers.

    Based on this, we can determine the attackers and which creatures will die.
    """

    # Creatures to attack with to force the optimal outcome
    attacking_creatures = [1, 2, 3, 4]

    # Creatures that will die in the resulting combat
    # (4) Wind Drake trades with (8) Skywinder Drake.
    # (5) River Bear (2/2) dies blocking (2) Centaur Courser (3/3).
    # (6) Grizzly Bears (2/2) dies blocking (3) Axebane Beast (2/5).
    # (7) Rusted Sentinel (3/4) dies blocking (1) Axegrinder Giant (6/4).
    dying_creatures = [4, 5, 6, 7, 8]

    # Sort the lists in increasing order as required
    attacking_creatures.sort()
    dying_creatures.sort()

    # Format the output string as per the example
    attackers_str = ", ".join([f"({c})" for c in attacking_creatures])
    dying_str = ", ".join([f"({c})" for c in dying_creatures])

    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()