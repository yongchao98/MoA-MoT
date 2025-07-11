import sys

def solve_magic_puzzle():
    """
    Calculates the maximum life Player B will lose in one turn.
    """
    # --- Step 1: Analyze Initial State & Pre-Combat Main Phase ---
    # Player A has 2 Island, 2 Swamp, 2 Mountain (6 mana total).
    # Player A casts Iron Apprentice for {1} mana.
    # The Replication Specialist ability triggers. Player A pays {2} mana to create a token copy.
    # Both the original and the token Iron Apprentice are 1/1 artifact creatures.
    # Mana spent: 3. Mana remaining: 3.
    
    print("Player A's optimal turn:")
    print("-" * 40)
    print("1. Pre-Combat Main Phase:")
    print("  - Player A casts Iron Apprentice (it's a 1/1 creature).")
    print("  - Player A uses Replication Specialist's ability to create a token copy, which is also a 1/1.")
    
    # --- Step 2: Beginning of Combat ---
    # Mukotai Soulripper's ability triggers. Player A sacrifices the original Iron Apprentice.
    # The sacrificed Iron Apprentice's death trigger puts its +1/+1 counter on Ironhoof Boar.
    
    ironhoof_boar_initial_power = 4
    ironhoof_boar_buffed_power = ironhoof_boar_initial_power + 1  # Becomes 5/3
    
    # Mukotai Soulripper gets +2/+2 and first strike from its own ability.
    mukotai_soulripper_initial_power = 2
    mukotai_soulripper_buffed_power = mukotai_soulripper_initial_power + 2 # Becomes 4/4
    
    print("\n2. Beginning of Combat:")
    print("  - Player A sacrifices the original Iron Apprentice to Mukotai Soulripper's ability.")
    print(f"  - The Apprentice's counter moves to Ironhoof Boar, making it a {ironhoof_boar_buffed_power}/3.")
    print(f"  - Mukotai Soulripper gets +2/+2, becoming a {mukotai_soulripper_buffed_power}/4 with menace and first strike.")
    
    # --- Step 3: Combat Phase (Declare Attackers) ---
    # Player A attacks with all available creatures.
    scrap_welder_power = 3
    replication_specialist_power = 2
    iron_apprentice_token_power = 1
    ninja_token_power = 1
    
    # Player B's optimal block to minimize life loss:
    # Player B has two 1/1 flying blockers. They can block non-flying creatures.
    # To prevent the most damage, they should block the highest-power non-trampling attackers.
    # These are Scrap Welder (3 power) and Replication Specialist (2 power).
    damage_prevented = scrap_welder_power + replication_specialist_power

    print("\n3. Combat Phase:")
    print("  - Player A attacks with all creatures.")
    print("  - To minimize life loss, Player B blocks Scrap Welder (preventing 3 damage) and Replication Specialist (preventing 2 damage).")
    print("\n  - The following damage gets through to Player B:")
    
    unblocked_creatures = {
        "Mukotai Soulripper": mukotai_soulripper_buffed_power,
        "Ironhoof Boar": ironhoof_boar_buffed_power,
        "Iron Apprentice Token": iron_apprentice_token_power,
        "Ninja Token #1": ninja_token_power,
        "Ninja Token #2": ninja_token_power,
    }

    combat_damage = 0
    for name, power in unblocked_creatures.items():
        print(f"    - {name}: {power} damage")
        combat_damage += power

    # --- Step 4: Post-Combat Main Phase ---
    # Player A uses the remaining mana to activate the Channel ability of Twinshot Sniper.
    channel_damage = 2
    print(f"\n4. Post-Combat Main Phase:")
    print(f"  - Player A uses the Channel ability of Twinshot Sniper to deal {channel_damage} direct damage.")
    
    # --- Step 5: Final Calculation ---
    total_life_lost = combat_damage + channel_damage
    
    print("-" * 40)
    print("\nTotal life lost by Player B this turn:")
    
    # Fulfilling the request to print each number in the final equation
    unblocked_powers = list(unblocked_creatures.values())
    equation = f"Total Damage = {unblocked_powers[0]} (Mukotai) + {unblocked_powers[1]} (Boar) + {unblocked_powers[2]} (Apprentice) + {unblocked_powers[3]} (Ninja) + {unblocked_powers[4]} (Ninja) + {channel_damage} (Channel)"
    result = f" = {total_life_lost}"
    
    print(equation + result)

solve_magic_puzzle()