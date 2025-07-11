def solve_mtg_puzzle():
    """
    Solves the MTG puzzle by determining the optimal attack and its outcome.
    """
    my_creatures = {
        1: {"name": "Axegrinder Giant", "power": 6, "toughness": 4, "can_attack": True},
        2: {"name": "Centaur Courser", "power": 3, "toughness": 3, "can_attack": True},
        3: {"name": "Axebane Beast", "power": 2, "toughness": 5, "can_attack": True},
        4: {"name": "Wind Drake", "power": 2, "toughness": 2, "can_attack": True, "flying": True},
    }

    opponent_creatures = {
        5: {"name": "River Bear", "power": 2, "toughness": 2, "can_block": True},
        6: {"name": "Grizzly Bears", "power": 2, "toughness": 2, "can_block": True},
        7: {"name": "Rusted Sentinel", "power": 3, "toughness": 4, "can_block": False}, # Defender
        8: {"name": "Skywinder Drake", "power": 2, "toughness": 2, "can_block": False}, # Defender
    }
    
    opponent_life = 3

    # Step 1: Find potential attackers and blockers
    potential_attackers = {num: stats for num, stats in my_creatures.items() if stats.get("can_attack")}
    potential_blockers = {num: stats for num, stats in opponent_creatures.items() if stats.get("can_block")}
    
    # Step 2: Determine unblockable damage
    # Wind Drake has flying, and the opponent has no flying blockers without defender.
    unblockable_damage = my_creatures[4]["power"]
    
    # Step 3: To deal lethal damage (3 total), we need at least 1 more damage from ground creatures.
    # The opponent has 2 blockers. To force damage through, we must attack with more than 2 ground creatures.
    # Therefore, we must attack with all creatures.
    attacking_creatures_nums = sorted(potential_attackers.keys())

    ground_attackers = sorted(
        [num for num in attacking_creatures_nums if not my_creatures[num].get("flying")],
        key=lambda num: my_creatures[num]["power"],
        reverse=True
    )
    
    # Step 4: Opponent makes optimal blocks to survive. They block the strongest attackers.
    dying_creatures_nums = []
    
    # The two strongest ground attackers are blocked.
    blocked_attackers = ground_attackers[:len(potential_blockers)]
    blocker_nums = sorted(potential_blockers.keys())

    # Simulate combat for the blocked creatures
    for i, attacker_num in enumerate(blocked_attackers):
        blocker_num = blocker_nums[i]
        
        attacker = my_creatures[attacker_num]
        blocker = opponent_creatures[blocker_num]
        
        # Check if blocker dies
        if attacker["power"] >= blocker["toughness"]:
            dying_creatures_nums.append(blocker_num)
            
        # Check if attacker dies
        if blocker["power"] >= attacker["toughness"]:
            dying_creatures_nums.append(attacker_num)
            
    dying_creatures_nums.sort()
    
    # Formatting the output string as requested
    # We must explicitly write each number in the final string
    attackers_str = f"({attacking_creatures_nums[0]}), ({attacking_creatures_nums[1]}), ({attacking_creatures_nums[2]}), ({attacking_creatures_nums[3]})"
    dying_str = f"({dying_creatures_nums[0]}), ({dying_creatures_nums[1]})"

    print(f"{attackers_str}; {dying_str}")

solve_mtg_puzzle()