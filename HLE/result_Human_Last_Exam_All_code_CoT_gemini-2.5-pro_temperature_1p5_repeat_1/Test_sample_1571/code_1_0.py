def calculate_max_damage():
    """
    Calculates the maximum life Player B will lose based on the described Magic: The Gathering board state.
    """

    # --- Initial State & Plan ---
    print("Player A's optimal plan involves three main phases: removing blockers, enhancing the attack force, and then attacking.")
    print("-" * 20)

    # --- Phase 1: Removing Blockers ---
    print("Phase 1: Remove Player B's two 1/1 flying blockers.")
    print("Action 1: Use Twinshot Sniper's Channel ability ({1}{R}) to deal 2 damage to one flyer.")
    print("Action 2: Cast March of Wretched Sorrow (X=1), exiling Junji to reduce the cost to {B}, to kill the second flyer.")
    print("Result: Player B now has no creatures to block with.")
    print("-" * 20)

    # --- Phase 2: Enhancing the Attack Force ---
    print("Phase 2: Cast Iron Apprentice. It enters as a 2/2 Artifact Creature.")
    print("-" * 20)

    # --- Phase 3: Combat and Damage Calculation ---
    print("Phase 3: Declare attackers and resolve triggers.")
    print("Player A attacks with Mukotai Soulripper, Scrap Welder, Ironhoof Boar, and two Ninja tokens.")
    
    # Calculate each creature's final power for the attack
    
    # Mukotai Soulripper: Base 4/4. Sacrifices Iron Apprentice for its attack trigger.
    mukotai_base_power = 4
    mukotai_trigger_bonus = 1
    mukotai_final_power = mukotai_base_power + mukotai_trigger_bonus
    print(f"Mukotai Soulripper (base power {mukotai_base_power}) sacrifices Iron Apprentice to get +{mukotai_trigger_bonus}/+{mukotai_trigger_bonus}, becoming a {mukotai_final_power}/{mukotai_final_power}.")
    
    # Scrap Welder: Base 3/3, enchanted to 4/4. Gets Iron Apprentice's counters.
    scrap_welder_base_power = 3
    clawing_torment_bonus = 1
    iron_apprentice_counters = 2
    scrap_welder_final_power = scrap_welder_base_power + clawing_torment_bonus + iron_apprentice_counters
    print(f"Scrap Welder (base power {scrap_welder_base_power}, +{clawing_torment_bonus}/+{clawing_torment_bonus} from enchantment) gets {iron_apprentice_counters} +1/+1 counters from the sacrificed Iron Apprentice, becoming a {scrap_welder_final_power}/{scrap_welder_final_power}.")

    # Ironhoof Boar: Base 4/3. Gets bonus on attack.
    ironhoof_boar_base_power = 4
    boar_attack_bonus = 2
    ironhoof_boar_final_power = ironhoof_boar_base_power + boar_attack_bonus
    print(f"Ironhoof Boar (base power {ironhoof_boar_base_power}) gets +{boar_attack_bonus}/+0 from its attack trigger, becoming a {ironhoof_boar_final_power}/3.")
    
    # Ninja Tokens
    ninja_1_power = 1
    ninja_2_power = 1
    print(f"The two unblockable Ninja tokens attack with {ninja_1_power} power each.")
    print("-" * 20)

    # --- Final Calculation ---
    total_damage = mukotai_final_power + scrap_welder_final_power + ironhoof_boar_final_power + ninja_1_power + ninja_2_power
    
    print("Since Player B has no blockers, all damage is dealt.")
    print("The final calculation for total life lost is:")
    # The prompt requests printing each number in the final equation.
    print(f"{mukotai_final_power} (from Mukotai) + {scrap_welder_final_power} (from Scrap Welder) + {ironhoof_boar_final_power} (from Ironhoof Boar) + {ninja_1_power} (from Ninja) + {ninja_2_power} (from Ninja) = {total_damage}")
    
    return total_damage

if __name__ == '__main__':
    max_life_loss = calculate_max_damage()
    print(f"\nThe maximum amount of life Player B will lose is {max_life_loss}.")