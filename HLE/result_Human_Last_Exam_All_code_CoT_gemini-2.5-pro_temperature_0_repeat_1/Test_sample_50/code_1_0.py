def calculate_max_damage():
    """
    Calculates the maximum possible damage based on the D&D Time Stop scenario.
    """
    spellcasting_modifier = 5

    # Delayed Blast Fireball (7th-level slot)
    # Cast on Turn 1, explodes on Turn 3.
    # We assume 1 round passes after the first, increasing damage by 1d6.
    dbf_base_dice = 12
    dbf_bonus_dice = 1 # Number of extra rounds of concentration
    dbf_total_dice = dbf_base_dice + dbf_bonus_dice
    dbf_damage = dbf_total_dice * 6
    
    # Scenario 1: Disintegrate (8th) and Spiritual Weapon (6th)
    # Disintegrate upcast to 8th level
    disintegrate_lvl_8_dice = 10 + (8 - 6) * 3
    disintegrate_lvl_8_damage = disintegrate_lvl_8_dice * 6 + 40
    
    # Spiritual Weapon upcast to 6th level
    # Base 1d8, +1d8 for every 2 levels above 2nd. 6th is 2 slots (4 levels) above.
    sw_lvl_6_dice = 1 + (6 - 2) // 2
    sw_lvl_6_damage = sw_lvl_6_dice * 8 + spellcasting_modifier
    
    total_damage_scenario_1 = dbf_damage + disintegrate_lvl_8_damage + sw_lvl_6_damage

    # Scenario 2: Disintegrate (6th) and Spiritual Weapon (8th)
    # Disintegrate at 6th level
    disintegrate_lvl_6_dice = 10
    disintegrate_lvl_6_damage = disintegrate_lvl_6_dice * 6 + 40

    # Spiritual Weapon upcast to 8th level
    # 8th is 3 slots (6 levels) above 2nd.
    sw_lvl_8_dice = 1 + (8 - 2) // 2
    sw_lvl_8_damage = sw_lvl_8_dice * 8 + spellcasting_modifier

    total_damage_scenario_2 = dbf_damage + disintegrate_lvl_6_damage + sw_lvl_8_damage

    # Determine the best scenario and print the result
    if total_damage_scenario_1 > total_damage_scenario_2:
        final_damage = total_damage_scenario_1
        print("The optimal strategy is to use the 8th-level slot for Disintegrate and the 6th-level slot for Spiritual Weapon.")
        print(f"Delayed Blast Fireball (7th level, {dbf_total_dice}d6) damage: {dbf_damage}")
        print(f"Disintegrate (8th level, {disintegrate_lvl_8_dice}d6 + 40) damage: {disintegrate_lvl_8_damage}")
        print(f"Spiritual Weapon (6th level, {sw_lvl_6_dice}d8 + 5) damage: {sw_lvl_6_damage}")
        print(f"\nTotal Damage = {dbf_damage} + {disintegrate_lvl_8_damage} + {sw_lvl_6_damage} = {final_damage}")
    else:
        final_damage = total_damage_scenario_2
        print("The optimal strategy is to use the 8th-level slot for Spiritual Weapon and the 6th-level slot for Disintegrate.")
        print(f"Delayed Blast Fireball (7th level, {dbf_total_dice}d6) damage: {dbf_damage}")
        print(f"Disintegrate (6th level, {disintegrate_lvl_6_dice}d6 + 40) damage: {disintegrate_lvl_6_damage}")
        print(f"Spiritual Weapon (8th level, {sw_lvl_8_dice}d8 + 5) damage: {sw_lvl_8_damage}")
        print(f"\nTotal Damage = {dbf_damage} + {disintegrate_lvl_6_damage} + {sw_lvl_8_damage} = {final_damage}")

calculate_max_damage()