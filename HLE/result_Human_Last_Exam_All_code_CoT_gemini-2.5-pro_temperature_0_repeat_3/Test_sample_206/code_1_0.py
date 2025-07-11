def solve_mtg_puzzle():
    """
    This function determines the solution to the MTG puzzle and prints the answer.
    """
    
    # Step 1: Define the creatures you should attack with for the win.
    # To win, you must force a lethal attacker through. By attacking with the Giant,
    # the Courser, and the Drake, you force the opponent into a complex blocking
    # situation. The key is that they are compelled to double-block the Giant,
    # which leaves the lethal Courser unblocked.
    attacking_creatures = [1, 2, 4]
    
    # Step 2: Define the creatures that will die in the resulting combat.
    # (1) Axegrinder Giant dies to the double block from (5) and (6).
    # (5) River Bear and (6) Grizzly Bears die blocking the Giant.
    # (4) Wind Drake trades with (8) Skywinder Drake.
    # (2) Centaur Courser is unblocked and survives.
    # (3) Axebane Beast did not attack and survives.
    # (7) Rusted Sentinel was not used to block and survives.
    dying_creatures = [1, 4, 5, 6, 8]
    
    # Step 3: Format the output as specified in the example.
    # The lists are already in increasing order.
    attackers_str = ", ".join(f"({c})" for c in attacking_creatures)
    dying_str = ", ".join(f"({c})" for c in dying_creatures)
    
    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()