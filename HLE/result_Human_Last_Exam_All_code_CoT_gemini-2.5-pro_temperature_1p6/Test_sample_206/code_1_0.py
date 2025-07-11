def solve_mtg_puzzle():
    """
    This function determines the optimal attack and resulting casualties for the MTG scenario.
    """

    # Step 1: Define the attacking creatures.
    # The optimal attack is with Axegrinder Giant (1), Centaur Courser (2), and Wind Drake (4).
    # This forces the opponent to block with all their available creatures to survive.
    attackers = [1, 2, 4]
    attackers.sort()

    # Step 2: Determine which creatures die if the opponent blocks optimally.
    # (4) Wind Drake (2/2) is blocked by (8) Skywinder Drake (2/2). Both die.
    # (2) Centaur Courser (3/3) is blocked by (6) Grizzly Bears (2/2). Grizzly Bears dies.
    # (1) Axegrinder Giant (6/4) is blocked by (7) Rusted Sentinel (3/4). Rusted Sentinel dies.
    creatures_that_die = [4, 6, 7, 8]
    creatures_that_die.sort()

    # Step 3: Format the output as specified in the example.
    # The format is (attacker1), (attacker2); (dead1), (dead2)
    attacker_str = ", ".join(f"({n})" for n in attackers)
    dead_str = ", ".join(f"({n})" for n in creatures_that_die)

    print(f"{attacker_str}; {dead_str}")

solve_mtg_puzzle()