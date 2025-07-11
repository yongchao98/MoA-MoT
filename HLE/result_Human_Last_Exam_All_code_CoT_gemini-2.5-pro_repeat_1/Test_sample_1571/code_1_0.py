def calculate_damage():
    """
    Calculates the maximum life Player B will lose based on the optimal plays from both players.
    """
    print("Player A's optimal turn plan to maximize damage:")
    
    # Step 1: Pre-combat actions
    print("\n--- Pre-Combat Main Phase ---")
    print("1. Player A casts Iron Apprentice (1 mana) and pays 2 mana for Replication Specialist's trigger, creating a 2/2 token copy.")
    print("2. Player A activates Twinshot Sniper's Channel ability (2 mana).")
    
    channel_damage = 2
    print(f"   - Twinshot Sniper deals {channel_damage} direct damage to Player B.")

    # Step 2: Combat phase
    print("\n--- Combat Phase ---")
    print("1. Player A crews Mukotai Soulripper with two 1/1 Ninja tokens.")
    print("2. Player A attacks with Mukotai Soulripper, Ironhoof Boar, Replication Specialist, Scrap Welder, and two Iron Apprentices.")
    print("3. On attack, Player A sacrifices a Ninja to give Mukotai Soulripper +1/+1 and menace. It is now a 5/6.")

    # Step 3: Player B's optimal blocks
    print("\n--- Declare Blockers Phase ---")
    print("Player B has two 1/1 flying blockers and must minimize damage.")
    print("Player A's attackers are: Mukotai (5/6, menace), Boar (4/3), Specialist (3/4), Welder (3/3), Apprentice (2/2), Apprentice (2/2).")
    print("Player B's optimal move is to block the two strongest non-menace creatures:")
    print(" - Blocker 1 on Ironhoof Boar (prevents 4 damage).")
    print(" - Blocker 2 on Replication Specialist (prevents 3 damage).")

    # Step 4: Damage Calculation
    print("\n--- Damage Calculation ---")
    print("The following attackers are unblocked and deal combat damage:")
    mukotai_damage = 5
    welder_damage = 3
    apprentice1_damage = 2
    apprentice2_damage = 2
    
    print(f" - Mukotai Soulripper deals {mukotai_damage} damage.")
    print(f" - Scrap Welder deals {welder_damage} damage.")
    print(f" - Iron Apprentice deals {apprentice1_damage} damage.")
    print(f" - The other Iron Apprentice deals {apprentice2_damage} damage.")

    total_combat_damage = mukotai_damage + welder_damage + apprentice1_damage + apprentice2_damage
    total_damage = channel_damage + total_combat_damage

    print("\n--- Total Life Lost ---")
    print("The total life Player B loses is the sum of the direct damage and the unblocked combat damage.")
    print(f"{channel_damage} (Channel) + {mukotai_damage} (Mukotai) + {welder_damage} (Welder) + {apprentice1_damage} (Apprentice) + {apprentice2_damage} (Apprentice) = {total_damage}")

calculate_damage()