def calculate_damage():
    """
    Calculates the maximum life Player B will lose in one turn.
    """

    # --- Pre-Combat Main Phase ---
    print("Player A's Optimal Pre-Combat Actions:")
    
    # Player A casts Twinshot Sniper. Its ETB (Enters the Battlefield) ability deals 2 damage.
    twinshot_sniper_etb_damage = 2
    print(f"- Player A casts Twinshot Sniper. Its ability triggers, dealing {twinshot_sniper_etb_damage} damage.")

    # Casting an artifact triggers Replication Specialist. Player A pays {2} to create a copy.
    # The copy's ETB also deals 2 damage.
    copied_sniper_etb_damage = 2
    print(f"- Replication Specialist creates a copy of Twinshot Sniper. The copy's ability triggers, dealing another {copied_sniper_etb_damage} damage.")

    pre_combat_damage = twinshot_sniper_etb_damage + copied_sniper_etb_damage
    print(f"\nTotal Pre-Combat Damage: {pre_combat_damage}\n")

    # --- Combat Phase ---
    print("Player A's Combat Phase:")
    # Player A crews Mukotai Soulripper and sacrifices a Ninja to make it a 4/3 with Menace.
    # Then attacks with all creatures.
    attackers = {
        "Mukotai Soulripper": 4,
        "Ironhoof Boar": 4,
        "Scrap Welder": 2,
        "Twinshot Sniper": 2,
        "Twinshot Sniper Token": 2,
        "Replication Specialist": 1
    }
    print("Player A attacks with: Mukotai Soulripper (4/3 Menace), Ironhoof Boar (4/3), Scrap Welder (2/2), Replication Specialist (1/3), and two Twinshot Snipers (2/2).")

    print("\nPlayer B's Optimal Blocking:")
    # Player B wants to minimize damage. They have two 1/1 blockers.
    # Blocking the 4/3 Menace creature requires both blockers. Damage prevented: 4. Damage taken from others: 4+2+2+2+1 = 11.
    # Blocking the two highest-power non-menace creatures (Ironhoof Boar and a 2/2) is better.
    # Damage prevented: 4 (from Boar) + 2 (from a 2/2) = 6.
    # Damage taken from others: 4 (Mukotai) + 2 (Sniper) + 2 (Sniper Token) + 1 (Specialist) = 9.
    # Player B chooses to block the Ironhoof Boar and Scrap Welder.
    blocked_creatures = ["Ironhoof Boar", "Scrap Welder"]
    print(f"Player B uses their two 1/1 blockers to block the {blocked_creatures[0]} and the {blocked_creatures[1]}.")
    
    unblocked_attackers = {name: power for name, power in attackers.items() if name not in blocked_creatures}
    
    print("\nUnblocked creatures dealing damage:")
    combat_damage_sources = []
    for name, power in unblocked_attackers.items():
        print(f"- {name} deals {power} damage.")
        combat_damage_sources.append(power)
        
    combat_damage = sum(combat_damage_sources)
    print(f"\nTotal Combat Damage: {combat_damage}\n")

    # --- Total Life Loss ---
    total_damage = pre_combat_damage + combat_damage
    
    # Construct the final equation string
    all_damage_sources = [twinshot_sniper_etb_damage, copied_sniper_etb_damage] + combat_damage_sources
    equation = " + ".join(map(str, all_damage_sources))

    print("--- Final Calculation ---")
    print("The total life Player B will lose is the sum of the pre-combat and combat damage.")
    print(f"Final Equation: {equation} = {total_damage}")

calculate_damage()