def solve_mtg_puzzle():
    """
    This function calculates the optimal attack and the resulting dead creatures for the MTG puzzle.
    """
    # 1. Define the creatures that will participate in the attack.
    # The Axegrinder Giant (1) is forced to attack.
    # Attacking with the Wind Drake (4) and Centaur Courser (2) is optimal
    # to overwhelm the opponent's limited number of blockers.
    attacking_creatures = [1, 2, 4]
    attacking_creatures.sort()

    # 2. Determine which creatures will die based on optimal blocking.
    # The Rusted Sentinel (7) cannot block. The opponent has 3 blockers for 3 attackers.
    # To survive, they must block all three.
    dead_creatures = []

    # Block 1: The flying (4) Wind Drake (2/2) is blocked by the (8) Skywinder Drake (2/2).
    # They trade, and both die.
    dead_creatures.append(4)
    dead_creatures.append(8)

    # Block 2 & 3: The ground attackers (1) Axegrinder Giant (4/3) and (2) Centaur Courser (3/3)
    # are blocked by (5) River Bear (2/2) and (6) Grizzly Bears (2/2).
    # Both attackers have more power than the blockers' toughness, so both blockers die.
    dead_creatures.append(5)
    dead_creatures.append(6)

    dead_creatures.sort()

    # 3. Format the output string as requested.
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
    dead_creatures_str = ", ".join(f"({c})" for c in dead_creatures)
    
    print(f"{attackers_str}; {dead_creatures_str}")

solve_mtg_puzzle()