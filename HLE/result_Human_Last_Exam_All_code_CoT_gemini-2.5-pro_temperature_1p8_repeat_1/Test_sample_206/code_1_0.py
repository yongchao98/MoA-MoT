def solve_mtg_puzzle():
    """
    This function calculates and prints the solution to the MTG puzzle.
    
    The plan is as follows:
    1.  Determine the necessary action: Win this turn, which requires dealing 3 damage.
    2.  Identify the overload principle: Attacking with three ground creatures against two blockers guarantees one will get through.
    3.  Select the attack: Attack with all creatures (1, 2, 3, 4) to maximize pressure. The unblocked (3) Beast deals 2 damage, forcing the opponent to 1 life. This makes the attacking (4) Wind Drake a lethal threat, forcing a block.
    4.  Determine optimal blocks by the opponent:
        - The opponent must block the strongest ground attackers, (1) and (2), leaving (3) to deal 2 damage.
        - To survive, they must block the (4) Wind Drake with their (8) Skywinder Drake.
    5.  Calculate combat results to see which creatures die.
        - (4) Wind Drake (2/2) trades with (8) Skywinder Drake (3/1). Both die.
        - (2) Centaur Courser (3/3) is blocked by (7) Rusted Sentinel (3/4). The Courser dies, the Sentinel survives.
        - (1) Axegrinder Giant (6/4) is blocked by (6) Grizzly Bears (2/2). The Giant survives, the Bears die.
    6.  List the attacking creatures and the creatures that die, sorted in increasing order.
    """
    
    # Creatures you should attack with
    attacking_creatures = [1, 2, 3, 4]
    
    # Creatures that will die after the attack
    # (2) dies to the 3/4 Sentinel
    # (4) dies to the 3/1 Drake
    # (6) dies to the 6/4 Giant
    # (8) dies to the 2/2 Drake
    dying_creatures = [2, 4, 6, 8]
    
    attack_str = ", ".join(f"({c})" for c in attacking_creatures)
    died_str = ", ".join(f"({c})" for c in dying_creatures)
    
    print(f"{attack_str}; {died_str}")

solve_mtg_puzzle()