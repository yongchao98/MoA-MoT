def calculate_damage():
    """
    This function calculates the maximum life loss for Player B by following the optimal play-by-play.
    """
    
    # Step 1: Pre-combat direct damage from Twinshot Sniper's Channel ability
    twinshot_damage = 2
    print(f"Player A starts by using Twinshot Sniper's Channel ability.")
    print(f"This deals {twinshot_damage} direct damage to Player B.")
    
    # Step 2: Set up for combat
    print("\nPlayer A then casts Junji, the Midnight Sky and crews Mukotai Soulripper.")
    print("When Mukotai Soulripper attacks, Player A sacrifices a tapped Ninja to make it a 3/4.")
    print("Player B must use their two 1/1 flying blockers on the 5/5 Junji because of Menace.")

    # Step 3: Calculate unblocked combat damage
    print("\nThe following attackers are unblocked and deal combat damage:")
    soulripper_damage = 3
    boar_damage = 4
    specialist_damage = 2
    ninja_damage = 2
    welder_damage = 1

    print(f"- Mukotai Soulripper (3/4 Flying, Menace) deals: {soulripper_damage} damage")
    print(f"- Ironhoof Boar (4/3 Trample) deals: {boar_damage} damage")
    print(f"- Replication Specialist (2/4) deals: {specialist_damage} damage")
    print(f"- Untapped Ninja (2/2) deals: {ninja_damage} damage")
    print(f"- Scrap Welder (1/3) deals: {welder_damage} damage")
    
    # Step 4: Sum all damage sources
    combat_damage = soulripper_damage + boar_damage + specialist_damage + ninja_damage + welder_damage
    total_damage = twinshot_damage + combat_damage

    # Final Output
    print("\nThe total life Player B loses is the sum of the initial damage and the combat damage.")
    print(f"Final calculation: {twinshot_damage} + {soulripper_damage} + {boar_damage} + {specialist_damage} + {ninja_damage} + {welder_damage} = {total_damage}")
    print(f"\nTotal life lost by Player B: {total_damage}")

calculate_damage()