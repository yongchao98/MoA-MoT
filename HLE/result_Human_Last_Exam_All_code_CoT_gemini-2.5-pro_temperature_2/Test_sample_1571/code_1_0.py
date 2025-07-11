def solve_magic_puzzle():
    """
    Calculates the maximum life Player B will lose this turn.
    """
    # Player A's starting board state
    # All are untapped, so we only need their power.
    replication_specialist_power = 2
    scrap_welder_base_power = 3
    ironhoof_boar_power = 4
    ninja_token_power = 1
    num_ninja_tokens = 2
    
    # Player B's enchantments on Player A's creatures
    clawing_torment_power_reduction = -1
    scrap_welder_power = scrap_welder_base_power + clawing_torment_power_reduction
    
    # --- Player A's Main Phase ---
    # Player A has 6 mana and plays two creatures from hand.
    # 1. Cast Iron Apprentice (enters as a 1/1)
    iron_apprentice_power = 1
    # 2. Cast Junji, the Midnight Sky (5/5 Flying, Menace)
    junji_power = 5

    print("Player A's Plan:")
    print("1. Cast Iron Apprentice (1/1) for 1 mana.")
    print("2. Cast Junji, the Midnight Sky (5/5) for 5 mana.")
    print("-" * 20)

    # --- Player A's Combat Phase ---
    # 1. Crew Mukotai Soulripper with one 1/1 Ninja token.
    #    The ninja is now tapped and cannot attack.
    mukotai_soulripper_base_power = 2
    
    # 2. Declare all possible attackers.
    attackers = {
        "Junji": junji_power,
        "Replication Specialist": replication_specialist_power,
        "Scrap Welder": scrap_welder_power,
        "Ironhoof Boar": ironhoof_boar_power,
        "Ninja Token": ninja_token_power,
        "Iron Apprentice": iron_apprentice_power,
        "Mukotai Soulripper": mukotai_soulripper_base_power,
    }
    print("Player A declares attackers:", list(attackers.keys()))

    # 3. Resolve attack triggers to maximize damage.
    # Player A sacrifices Iron Apprentice to Mukotai Soulripper's ability.
    soulripper_sacrifice_bonus = 2
    attackers["Mukotai Soulripper"] += soulripper_sacrifice_bonus
    attackers.pop("Iron Apprentice") # Sacrificed creature is removed from combat
    print("Player A sacrifices Iron Apprentice to give Mukotai Soulripper +2/+0.")
    
    # Iron Apprentice's death trigger puts its +1/+1 counter on an artifact creature.
    # To maximize damage, the counter should go on an unblocked creature.
    # Player A places the counter on the now 4/2 Soulripper, making it a 5/3 with Menace.
    iron_apprentice_counter = 1
    attackers["Mukotai Soulripper"] += iron_apprentice_counter
    print("Iron Apprentice's +1/+1 counter is put on Mukotai Soulripper.")
    print("-" * 20)

    # --- Player B's Blocking Phase ---
    # Player B has 2x 1/1 Mothrider Patrols. They must block to minimize life loss.
    # Junji has flying, so they cannot block it. It deals 5 damage.
    unblocked_attackers = {}
    unblocked_attackers["Junji (Flying)"] = attackers.pop("Junji")
    
    # The remaining attackers are on the ground. Player B must block the two with the highest power.
    # Mukotai Soulripper has Menace, so it cannot be blocked by a single creature.
    # The highest-power non-Menace threats are Ironhoof Boar (4) and Replication Specialist (2).
    # Player B blocks them.
    blocked_attackers = {}
    blocked_attackers["Ironhoof Boar"] = attackers.pop("Ironhoof Boar")
    blocked_attackers["Replication Specialist"] = attackers.pop("Replication Specialist")

    print("Player B's optimal block to minimize life loss:")
    print(" - Junji is unblockable (Flying).")
    print(" - Player B's two 1/1 creatures block the two strongest non-menace attackers:")
    print(f"   - One blocks the {list(blocked_attackers.keys())[0]} ({blocked_attackers['Ironhoof Boar']}).")
    print(f"   - The other blocks the {list(blocked_attackers.keys())[1]} ({blocked_attackers['Replication Specialist']}).")
    print("-" * 20)

    # The rest of the attackers are unblocked.
    unblocked_attackers.update(attackers)
    
    print("The following attackers are unblocked:")
    total_damage = 0
    damage_sources = []
    for name, power in unblocked_attackers.items():
        print(f" - {name}: {power} damage")
        total_damage += power
        damage_sources.append(str(power))

    # --- Final Calculation ---
    print("\nFinal Life Loss Calculation:")
    print(" + ".join(damage_sources) + f" = {total_damage}")
    
    return total_damage

final_answer = solve_magic_puzzle()
# print(f"\n<<<The maximum life Player B will lose this turn is {final_answer}.>>>")