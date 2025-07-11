def calculate_life_loss():
    """
    Calculates the maximum life Player B will lose in one turn.
    """
    total_life_loss = 0
    damage_sources = []

    # Step 1: Pre-combat actions to maximize damage
    # Player A has 2 Island, 2 Swamp, 2 Mountain (2U, 2B, 2R)
    # Player A uses the Channel ability of Twinshot Sniper for {1}{R}.
    channel_damage = 2
    total_life_loss += channel_damage
    damage_sources.append(str(channel_damage))
    print(f"Player A starts by paying {{1}}{{R}} to use the Channel ability of Twinshot Sniper, dealing {channel_damage} damage to Player B.")
    
    # After Channel, A has 2U, 2B, 1R mana available.
    # Player A casts Iron Apprentice for {1} (using 1 Island).
    # This triggers Replication Specialist. A pays {2} (using 1 Island, 1 Mountain) to create a copy.
    # Player A now has two 1/1 Iron Apprentice creatures.
    # Remaining mana: 2B
    print("Player A casts Iron Apprentice for {1}, then pays {2} for Replication Specialist's triggered ability to create a token copy. They now control two 1/1 Iron Apprentices.")

    # Step 2: Crewing the vehicle
    # Player A taps the 2/4 Replication Specialist to crew the 4/4 Mukotai Soulripper.
    print("Player A crews the Mukotai Soulripper by tapping the Replication Specialist, turning it into a 4/4 creature with Menace.")

    # Step 3: Combat Phase
    # Player B has two 1/1 Mothrider Patrols to block.
    # Because Mukotai Soulripper has Menace, it must be blocked by two or more creatures.
    # To minimize damage, Player B must use both 1/1 creatures to block the 4/4 Mukotai Soulripper.
    # This prevents 4 damage, but leaves all other attackers unblocked.
    print("\nIn combat, Player B must use both of their 1/1 Mothrider Patrols to block the 4/4 Mukotai Soulripper due to Menace.")
    print("This means all other attackers are unblocked.")
    
    # Unblocked attackers and their damage
    unblocked_attackers = {
        "Scrap Welder": 3,
        "Ironhoof Boar (Haste)": 4,
        "Ninja token 1": 1,
        "Ninja token 2": 1,
        "Iron Apprentice 1": 1,
        "Iron Apprentice 2 (token)": 1,
    }

    print("\nDamage from unblocked creatures:")
    combat_damage = 0
    for creature, power in unblocked_attackers.items():
        print(f"- {creature}: {power} damage")
        combat_damage += power
        damage_sources.append(str(power))
    
    total_life_loss += combat_damage

    # Step 4: Final Calculation
    print("\n------------------------------------")
    print("Total life lost by Player B this turn:")
    equation = " + ".join(damage_sources)
    print(f"{equation} = {total_life_loss}")

calculate_life_loss()