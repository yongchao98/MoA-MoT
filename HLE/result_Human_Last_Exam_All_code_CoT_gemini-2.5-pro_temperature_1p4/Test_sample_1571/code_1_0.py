def calculate_max_damage():
    """
    Calculates the maximum life Player B will lose in one turn.
    """
    # 1. Pre-combat phase
    # Player A uses the Channel ability of Twinshot Sniper.
    sniper_channel_damage = 2
    direct_damage = sniper_channel_damage

    # 2. Combat phase
    # Player A's creatures attack.
    # Mukotai Soulripper's power is calculated after its triggers resolve.
    mukotai_base_power = 2
    sacrificed_apprentice_counters = 2
    mukotai_trigger_buff = 2
    mukotai_final_power = mukotai_base_power + sacrificed_apprentice_counters + mukotai_trigger_buff

    # Power of all attacking creatures
    # Note: One Iron Apprentice is sacrificed, the two Ninjas are tapped to crew.
    # Scrap Welder is a 0/1 due to Clawing Torment.
    attackers_power = {
        "Mukotai Soulripper (menace)": mukotai_final_power,
        "Ironhoof Boar": 4,
        "Replication Specialist": 2,
        "Iron Apprentice": 2,
        "Scrap Welder": 0
    }
    total_attacking_power = sum(attackers_power.values())

    # 3. Blocking phase
    # Player B optimally blocks to minimize damage. They have two 1/1 blockers.
    # They can either double-block the 6-power menace creature or single-block the two
    # strongest non-menace creatures (4-power and 2-power). Both prevent 6 damage.
    damage_prevented_by_blocking_menace = attackers_power["Mukotai Soulripper (menace)"]
    
    non_menace_powers = sorted([p for name, p in attackers_power.items() if "menace" not in name], reverse=True)
    damage_prevented_by_blocking_others = non_menace_powers[0] + non_menace_powers[1]
    
    max_damage_prevented = max(damage_prevented_by_blocking_menace, damage_prevented_by_blocking_others)
    
    combat_damage = total_attacking_power - max_damage_prevented

    # 4. Final calculation
    total_damage = direct_damage + combat_damage

    # Print the step-by-step reasoning
    print("To maximize damage, Player A takes the following steps:")
    print(f"\n1. Pre-combat, Player A deals {sniper_channel_damage} direct damage using Twinshot Sniper's Channel ability.")
    
    print("\n2. In combat, after triggers, the total power of attacking creatures is calculated:")
    power_breakdown = " + ".join(map(str, attackers_power.values()))
    print(f"   Attacking Power = {power_breakdown} = {total_attacking_power}")

    print("\n3. Player B blocks optimally to prevent the most damage.")
    print(f"   Maximum damage Player B can prevent = {max_damage_prevented}")
    
    print(f"\n4. Total combat damage is {total_attacking_power} (total power) - {max_damage_prevented} (blocked power) = {combat_damage}.")

    print("\nFinal Result:")
    print("The total life Player B will lose is:")
    print(f"{sniper_channel_damage} + {combat_damage} = {total_damage}")


calculate_max_damage()
<<<10>>>