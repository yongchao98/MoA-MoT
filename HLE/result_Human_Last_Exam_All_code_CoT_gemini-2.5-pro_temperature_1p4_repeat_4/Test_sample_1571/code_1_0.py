def solve():
    """
    Calculates the maximum life Player B will lose based on the optimal sequence of plays.
    """
    # Damage from Twinshot Sniper's Channel ability
    channel_damage = 2

    # Player B uses their two 1/1 blockers to block the 5/5 Junji to prevent the most damage.
    # The following creatures are unblocked:

    # Damage from the unblocked 4/4 Mukotai Soulripper
    mukotai_damage = 4

    # Damage from the unblocked 3/3 Scrap Welder (which has menace)
    scrap_welder_damage = 3

    # Damage from the unblocked 4/4 Ironhoof Boar
    ironhoof_boar_damage = 4

    # Damage from the unblocked 2/4 Replication Specialist
    replication_specialist_damage = 2
    
    total_damage = channel_damage + mukotai_damage + scrap_welder_damage + ironhoof_boar_damage + replication_specialist_damage

    print("Player A casts Junji, the Midnight Sky, then uses Twinshot Sniper's Channel ability for 2 damage.")
    print("Player A attacks with Junji (5/5 menace), Mukotai Soulripper (4/4), Scrap Welder (3/3 menace), Ironhoof Boar (4/4), and Replication Specialist (2/4).")
    print("To minimize life loss, Player B uses their two blockers to block the 5/5 Junji.")
    print("The remaining unblocked creatures deal combat damage.")
    print(f"Total life lost = (Twinshot Damage) + (Mukotai Damage) + (Scrap Welder Damage) + (Ironhoof Boar Damage) + (Replication Specialist Damage)")
    print(f"Calculation: {channel_damage} + {mukotai_damage} + {scrap_welder_damage} + {ironhoof_boar_damage} + {replication_specialist_damage} = {total_damage}")

solve()