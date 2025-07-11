def calculate_max_damage():
    """
    Calculates the maximum life Player B will lose during Player A's turn by
    determining the optimal sequence of plays for Player A and the optimal
    blocking strategy for Player B.
    """
    
    # 1. Pre-combat Phase Actions & Damage
    twinshot_channel_damage = 2
    
    print("Player A's Turn Breakdown:")
    print("----------------------------")
    print("\nI. Pre-Combat Main Phase:")
    print("1. Player A activates the Channel ability of Twinshot Sniper. This deals 2 damage to Player B.")
    print("2. Player A casts Iron Apprentice, triggering Replication Specialist. Player A pays the cost to create a token copy. Player A now controls two 1/1 Iron Apprentice creatures.")
    
    # 2. Combat Phase Actions
    print("\nII. Combat Phase:")
    print("1. Player A taps one 1/1 Iron Apprentice to Crew the 4/4 Mukotai Soulripper.")
    print("2. Player A attacks with Mukotai Soulripper, Ironhoof Boar, 2 Ninja tokens, Replication Specialist, Scrap Welder, and the other 1/1 Iron Apprentice.")
    print("3. Mukotai Soulripper's attack trigger resolves. Player A sacrifices the tapped Iron Apprentice that was used for crew.")
    print("   - This gives Mukotai Soulripper 'menace', meaning it can only be blocked by two or more creatures.")
    print("   - The sacrificed Iron Apprentice's ability puts a +1/+1 counter on Ironhoof Boar, making it a 5/4 with trample.")

    # 3. Player B's Optimal Blocks and Damage Calculation
    print("\nIII. Damage Calculation:")
    print("Player B has two 1/1 creatures to block. Their optimal strategy to minimize damage is to let the menace creature through and block the two non-trampling creatures with the highest power.")
    print("Player B blocks the 2/3 Replication Specialist and one 2/2 Ninja token. These blocked creatures deal no damage.")
    
    print("\nThe unblocked creatures deal the following combat damage:")
    mukotai_damage = 4
    boar_damage = 5
    welder_damage = 1
    ninja_damage = 2
    apprentice_damage = 1
    
    print(f"- Unblocked Mukotai Soulripper (menace): {mukotai_damage} damage")
    print(f"- Unblocked Ironhoof Boar (5/4 trample): {boar_damage} damage")
    print(f"- Unblocked Scrap Welder: {welder_damage} damage")
    print(f"- Unblocked Ninja Token: {ninja_damage} damage")
    print(f"- Unblocked Iron Apprentice: {apprentice_damage} damage")

    # 4. Final Total
    total_damage = twinshot_channel_damage + mukotai_damage + boar_damage + welder_damage + ninja_damage + apprentice_damage

    print("\nIV. Total Life Lost:")
    print("The total life Player B loses is the sum of the direct damage and combat damage.")
    print("\nFinal Equation:")
    print(f"{twinshot_channel_damage} (from Twinshot Sniper) + {mukotai_damage} (from Mukotai Soulripper) + {boar_damage} (from Ironhoof Boar) + {welder_damage} (from Scrap Welder) + {ninja_damage} (from Ninja) + {apprentice_damage} (from Iron Apprentice) = {total_damage}")


calculate_max_damage()