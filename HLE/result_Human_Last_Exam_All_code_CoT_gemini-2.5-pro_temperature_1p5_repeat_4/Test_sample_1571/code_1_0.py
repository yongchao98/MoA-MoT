def solve_magic_puzzle():
    """
    Calculates the maximum amount of life Player B will lose in one turn.
    """

    # Step 1: Pre-combat direct damage
    # Player A has 6 mana (UUBBRR). They spend {1}{R} to channel Twinshot Sniper.
    direct_damage = 2
    print(f"First, Player A uses the Channel ability of Twinshot Sniper, dealing {direct_damage} damage to Player B.")
    print("-" * 20)

    # Step 2: Set up and declare attackers
    # Player A spends {1} to cast Iron Apprentice and {2} for Replication Specialist's trigger.
    # They crew Mukotai Soulripper with two 1/1 Ninjas.
    # On attack, A sacrifices an Iron Apprentice token to boost Mukotai Soulripper.
    attackers = {
        "Mukotai Soulripper": {"power": 4, "is_blocked": False, "note": "(4/3 with Menace after sacrifice)"},
        "Ironhoof Boar": {"power": 4, "is_blocked": False, "note": "(4/3 with Trample)"},
        "Scrap Welder": {"power": 3, "is_blocked": False, "note": "(3/3)"},
        "Replication Specialist": {"power": 2, "is_blocked": False, "note": "(2/4)"},
        "Iron Apprentice": {"power": 1, "is_blocked": False, "note": "(1/1)"},
    }
    print("Player A declares the following attackers:")
    for name, data in attackers.items():
        print(f"- {name}: {data['power']} power {data['note']}")
    print("-" * 20)

    # Step 3: Player B's optimal blocking strategy
    # Player B has two 1/1 blockers. To minimize damage, they block the highest power
    # non-trample/non-menace creatures they can. These are Scrap Welder (3) and Replication Specialist (2).
    # Blocking Mukotai would take both blockers and let 10 damage through. This way only lets 9 through.
    attackers["Scrap Welder"]["is_blocked"] = True
    attackers["Replication Specialist"]["is_blocked"] = True
    print("Player B has two 1/1 blockers and must minimize life loss.")
    print("Player B's optimal strategy is to block the Scrap Welder and the Replication Specialist.")
    print("-" * 20)

    # Step 4: Calculate combat damage
    combat_damage = 0
    combat_damage_calculation = []
    print("The following attackers deal combat damage:")
    for name, data in attackers.items():
        if not data["is_blocked"]:
            damage = data['power']
            combat_damage += damage
            combat_damage_calculation.append(str(damage))
            print(f"- {name} deals {damage} damage.")

    print("\nTotal combat damage calculation:")
    combat_damage_equation = " + ".join(combat_damage_calculation)
    print(f"{combat_damage_equation} = {combat_damage}")
    print("-" * 20)
    
    # Step 5: Final Calculation
    total_life_loss = direct_damage + combat_damage
    print("Total life loss for Player B is the direct damage plus the combat damage.")
    print(f"Final Equation: {direct_damage} + {combat_damage} = {total_life_loss}")

solve_magic_puzzle()
<<<11>>>