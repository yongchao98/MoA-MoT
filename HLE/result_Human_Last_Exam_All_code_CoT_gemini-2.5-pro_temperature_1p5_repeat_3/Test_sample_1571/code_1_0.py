def calculate_max_damage():
    """
    Calculates the maximum life Player B can lose in a single turn based on the provided game state.
    """
    # Player A's creatures' initial power
    initial_mukotai_power = 3
    initial_scrap_welder_power = 3
    initial_ironhoof_boar_power = 4
    initial_replication_specialist_power = 2

    # Player B's life loss will come from direct damage and combat damage.

    print("To maximize the life Player B loses, Player A should execute the following plan:")
    print("-" * 70)

    # 1. Pre-combat Direct Damage
    print("Step 1: Pre-Combat Main Phase Actions")
    print("Player A spends {1}{R} to use the Channel ability of Twinshot Sniper.")
    direct_damage = 2
    print(f"This deals {direct_damage} damage directly to Player B. Twinshot Sniper goes to the graveyard.")
    print("Player A then casts Iron Apprentice for {1} mana. It enters as a 1/1 artifact creature.")
    print("-" * 70)

    # 2. Combat Phase
    print("Step 2: Combat Phase")
    print("Player A crews Mukotai Soulripper with their two 1/1 Ninja tokens.")
    print("Player A attacks with all available creatures: Mukotai Soulripper, Scrap Welder, Ironhoof Boar, Replication Specialist, and Iron Apprentice.")
    print("\nAttack triggers resolve:")
    
    # Mukotai trigger
    mukotai_bonus = 2
    final_mukotai_power = initial_mukotai_power + mukotai_bonus
    print(f"- Mukotai Soulripper's trigger: Player A sacrifices the Iron Apprentice to give it +{mukotai_bonus}/+0. Its power becomes {final_mukotai_power}.")
    
    # Iron Apprentice trigger
    scrap_welder_bonus = 1
    final_scrap_welder_power = initial_scrap_welder_power + scrap_welder_bonus
    print(f"- Iron Apprentice's death trigger: Player A puts its +1/+1 counter on Scrap Welder. Its power becomes {final_scrap_welder_power}.")
    print("-" * 70)

    # 3. Damage Calculation
    print("Step 3: Damage Calculation")
    print("Player B has two 1/1 creatures with Flying, but none of Player A's attackers have flying.")
    print("Therefore, Player B cannot block any of the attacking creatures.")
    
    combat_damage = final_mukotai_power + final_scrap_welder_power + initial_ironhoof_boar_power + initial_replication_specialist_power
    print(f"Unblocked combat damage is the sum of the attackers' power: {final_mukotai_power} + {final_scrap_welder_power} + {initial_ironhoof_boar_power} + {initial_replication_specialist_power} = {combat_damage}.")
    print("-" * 70)

    # 4. Total Life Loss
    print("Step 4: Final Life Loss Calculation")
    total_life_loss = direct_damage + combat_damage
    print("The total life Player B loses is the sum of the direct damage and combat damage.")
    print(f"Final Equation: {direct_damage} (from Twinshot Sniper) + {final_mukotai_power} (from Mukotai Soulripper) + {final_scrap_welder_power} (from Scrap Welder) + {initial_ironhoof_boar_power} (from Ironhoof Boar) + {initial_replication_specialist_power} (from Replication Specialist) = {total_life_loss}")

calculate_max_damage()