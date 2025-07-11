def solve_magic_puzzle():
    """
    Calculates the maximum life Player B will lose in the given Magic: The Gathering scenario.
    """

    # Step 1: Direct damage from pre-combat actions
    twinshot_sniper_channel_damage = 2
    print(f"Player A uses Twinshot Sniper's Channel ability, dealing {twinshot_sniper_channel_damage} direct damage to Player B.")

    # Step 2: Calculate combat damage
    # Player A's attackers after pre-combat setup. Player B has no legal blockers.
    replication_specialist_damage = 2
    # Mukotai is 4/5, crewed, attacks, sacs a ninja, becomes 6/5
    mukotai_soulripper_damage = 6
    scrap_welder_damage = 3
    ironhoof_boar_damage = 4
    # Iron Apprentice is cast and copied, both enter as 2/3s
    iron_apprentice_damage = 2
    iron_apprentice_token_damage = 2

    print("\nPlayer A attacks. Player B's flying creatures cannot block any of the ground attackers.")
    print(f"Combat damage from Replication Specialist: {replication_specialist_damage}")
    print(f"Combat damage from Mukotai Soulripper (after sacrificing a Ninja): {mukotai_soulripper_damage}")
    print(f"Combat damage from Scrap Welder: {scrap_welder_damage}")
    print(f"Combat damage from Ironhoof Boar: {ironhoof_boar_damage}")
    print(f"Combat damage from Iron Apprentice: {iron_apprentice_damage}")
    print(f"Combat damage from the Iron Apprentice token: {iron_apprentice_token_damage}")

    # Step 3: Sum all damage
    total_combat_damage = (replication_specialist_damage +
                           mukotai_soulripper_damage +
                           scrap_welder_damage +
                           ironhoof_boar_damage +
                           iron_apprentice_damage +
                           iron_apprentice_token_damage)

    total_life_loss = twinshot_sniper_channel_damage + total_combat_damage

    print(f"\nTotal life loss for Player B is the sum of all damage sources:")
    print(f"{twinshot_sniper_channel_damage} (Channel) + {replication_specialist_damage} + {mukotai_soulripper_damage} + {scrap_welder_damage} + {ironhoof_boar_damage} + {iron_apprentice_damage} + {iron_apprentice_token_damage} = {total_life_loss}")

solve_magic_puzzle()