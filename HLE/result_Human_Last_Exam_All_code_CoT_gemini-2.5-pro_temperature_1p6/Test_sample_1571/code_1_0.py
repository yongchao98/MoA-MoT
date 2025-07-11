def calculate_max_damage():
    """
    Calculates the maximum life Player B will lose in one turn.
    """
    life_lost = 0
    
    # Step 1: Pre-Combat Main Phase
    # Player A has 2 Island, 2 Swamp, and 2 Mountain, for a total of 6 mana.
    # The optimal play is to cast "March of Wretched Sorrow".
    # The cost is {X}{B}{B}. With 6 mana available, Player A can set X to 4,
    # paying {4}{B}{B} by tapping all their lands.
    # This spell causes Player B to lose X life.
    march_damage = 4
    life_lost += march_damage
    print(f"Step 1: Player A casts March of Wretched Sorrow for X={march_damage}.")
    print(f"Player B loses {march_damage} life. Total life lost: {life_lost}")
    print("---")
    
    # Step 2: Combat Phase - Attack Triggers
    # Player A crews Mukotai Soulripper (Crew 4) by tapping Ironhoof Boar (4 power).
    # Player A attacks with Mukotai Soulripper (4/5), Replication Specialist (2/4),
    # the two Ninja tokens (2/2 each), and Scrap Welder (currently a 2/2 due to Clawing Torment's -1/-1).
    # Clawing Torment has an ability: "Whenever this creature attacks, defending player loses 1 life."
    # When Scrap Welder attacks, this triggers.
    clawing_torment_trigger_damage = 1
    life_lost += clawing_torment_trigger_damage
    print(f"Step 2: Player A attacks with Scrap Welder, triggering Clawing Torment.")
    print(f"Player B loses {clawing_torment_trigger_damage} life. Total life lost: {march_damage} + {clawing_torment_trigger_damage} = {life_lost}")
    print("---")

    # Step 3: Combat Phase - Damage Calculation
    # Player B has two 2/1 Mothrider Patrols to block. To minimize damage, they must block the highest power attackers.
    # Attacking creatures' powers: Mukotai Soulripper (4), Replication Specialist (2), Scrap Welder (2), Ninja (2), Ninja (2).
    # Blocker 1 blocks the 4-power Mukotai Soulripper, preventing 4 damage.
    # Blocker 2 blocks one of the 2-power creatures (e.g., Scrap Welder), preventing 2 damage.
    # The remaining unblocked creatures deal their damage.
    unblocked_creature_1_damage = 2  # Replication Specialist
    unblocked_creature_2_damage = 2  # Ninja token
    unblocked_creature_3_damage = 2  # Ninja token
    combat_damage = unblocked_creature_1_damage + unblocked_creature_2_damage + unblocked_creature_3_damage
    life_lost += combat_damage
    
    print("Step 3: Player B blocks Mukotai Soulripper (4 damage) and Scrap Welder (2 damage).")
    print(f"Unblocked attackers (Replication Specialist, two Ninjas) deal combat damage.")
    print(f"Player B loses {unblocked_creature_1_damage} + {unblocked_creature_2_damage} + {unblocked_creature_3_damage} = {combat_damage} combat damage.")
    print("---")

    # Final Calculation
    print("Final Calculation:")
    print(f"Total Life Lost = (March of Wretched Sorrow) + (Clawing Torment trigger) + (Combat Damage)")
    print(f"Total Life Lost = {march_damage} + {clawing_torment_trigger_damage} + {combat_damage}")
    print(f"Final Answer: {life_lost}")

calculate_max_damage()