def solve_mtg_puzzle():
    """
    Solves the MTG puzzle by determining the optimal attack and resulting casualties.
    """
    
    # Define player and opponent creatures with their number, name, power, and toughness.
    # Your creatures
    creature_1 = {"id": 1, "name": "Axegrinder Giant", "power": 6, "toughness": 4}
    creature_2 = {"id": 2, "name": "Centaur Courser", "power": 3, "toughness": 3}
    
    # Opponent's creatures that can block
    creature_6 = {"id": 6, "name": "Grizzly Bears", "power": 2, "toughness": 2}
    creature_7 = {"id": 7, "name": "Rusted Sentinel", "power": 3, "toughness": 4}

    # The optimal attack is with creatures (1) and (2).
    attacking_creatures_ids = [creature_1["id"], creature_2["id"]]
    attacking_creatures_ids.sort()
    
    # Determine which creatures die in combat.
    # The opponent must block both lethal attackers to survive.
    # Optimal blocking is (7) on (1) and (6) on (2).
    
    dead_creatures_ids = []
    
    # Combat: (1) Axegrinder Giant vs (7) Rusted Sentinel
    # Giant's power (6) > Sentinel's toughness (4), so Sentinel dies.
    if creature_1["power"] >= creature_7["toughness"]:
        dead_creatures_ids.append(creature_7["id"])
    
    # Combat: (2) Centaur Courser vs (6) Grizzly Bears
    # Courser's power (3) > Bears' toughness (2), so Bears die.
    if creature_2["power"] >= creature_6["toughness"]:
        dead_creatures_ids.append(creature_6["id"])
        
    dead_creatures_ids.sort()

    # Format the output as requested.
    attackers_str = ", ".join(f"({id})" for id in attacking_creatures_ids)
    dead_str = ", ".join(f"({id})" for id in dead_creatures_ids)
    
    print(f"{attackers_str}; {dead_str}")

solve_mtg_puzzle()