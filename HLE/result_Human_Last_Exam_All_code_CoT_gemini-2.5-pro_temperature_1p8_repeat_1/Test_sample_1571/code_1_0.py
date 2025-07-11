def calculate_max_damage():
    """
    Calculates the maximum life Player B will lose in one turn.
    """
    # -- Initial State --
    player_b_life_loss = 0
    print("Player A's optimal line of play:")
    print("--------------------------------\n")

    # -- Pre-Combat Main Phase --
    print("1. Cast Iron Apprentice for {1} mana.")
    print("2. Replication Specialist triggers. Player A pays {2} to create a token copy of Iron Apprentice.")
    
    # Player A channels Twinshot Sniper for direct damage.
    direct_damage = 2
    player_b_life_loss += direct_damage
    print(f"3. Channel Twinshot Sniper for {{1}}{{R}}, dealing {direct_damage} damage to Player B.")

    # -- Combat Phase --
    # Attacking creatures' power after setup and triggers:
    # - Mukotai Soulripper: Base 4 power when crewed. On attack, a 1/1 Iron Apprentice is sacrificed.
    #   Mukotai gains +1 power. It becomes 5 power.
    # - Replication Specialist: Base 2 power. The sacrificed Iron Apprentice's death trigger gives it a +1/+1 counter.
    #   It becomes 3 power.
    # - Scrap Welder: Base 3 power, +1 from Clawing Torment. It has 4 power.
    # - Ironhoof Boar: 4 power.
    # - Remaining Ninja Token: 1 power.
    # - Remaining Iron Apprentice Token: 1 power.

    unblocked_attackers = {
        "Mukotai Soulripper": 5,
        "Replication Specialist": 3,
        "Ninja Token": 1,
        "Iron Apprentice Token": 1
    }

    print("\n4. Enter Combat.")
    print("   - Player A crews Mukotai Soulripper.")
    print("   - Player A attacks with all creatures.")
    print("   - On attack, Player A sacrifices an Iron Apprentice to Mukotai Soulripper.")
    print("     - Mukotai becomes a 5/4 with menace.")
    print("     - Replication Specialist becomes a 3/5 from the Apprentice's death trigger.")
    
    print("\n5. Player B Declares Blocks.")
    print("   - To minimize damage, Player B must block Scrap Welder (4 damage + 4 life loss) and Ironhoof Boar (4 damage).")
    print("   - All other attackers are unblocked.\n")

    # -- Damage Calculation --
    print("6. Calculating Total Life Loss.")
    print("   - Direct Damage (Twinshot Sniper): 2")
    combat_damage = 0
    combat_damage_sources = []
    print("   - Unblocked Combat Damage:")
    for name, power in unblocked_attackers.items():
        print(f"     - {name}: {power} damage")
        combat_damage += power
        combat_damage_sources.append(power)

    # Sort combat damage sources for the final equation printout
    combat_damage_sources.sort(reverse=True)

    # Sum all damage for the final result
    total_life_loss = direct_damage + combat_damage
    
    final_equation_parts = [str(direct_damage)] + [str(p) for p in combat_damage_sources]
    final_equation = " + ".join(final_equation_parts)

    print("\n--- Final Result ---")
    print(f"The total amount of life Player B will lose is {final_equation} = {total_life_loss}.")

calculate_max_damage()
<<<12>>>