import math

def calculate_homm3_damage():
    """
    Calculates the damage dealt in a specific Heroes of Might and Magic III scenario.
    This function breaks down the calculation step-by-step and prints the final equation.
    """
    # Step 1: Define base stats and hero skills
    num_attackers = 3
    attacker_base_damage = 50
    attacker_base_attack = 30
    defender_base_defense = 21

    hero_attack_skill = 19
    hero_defense_skill = 1

    # Step 2: Define modifiers from spells and special abilities
    # Devils are in defensive stance: +20% to their Defense skill
    defend_stance_multiplier = 1.20
    # Archangels deal +50% damage against Devils
    racial_bonus_multiplier = 1.50
    # Devils are affected by Protection from Water (cast with Expert Water Magic): 30% damage reduction
    spell_reduction_multiplier = 1.0 - 0.30

    # Step 3: Calculate the attacker's total Attack skill
    # Total Attack = Creature's Base Attack + Hero's Attack Skill
    total_attack = attacker_base_attack + hero_attack_skill

    # Step 4: Calculate the defender's total Defense skill
    # The bonus from the "Defend" command is applied after adding the hero's Defense skill.
    # Total Defense = (Creature's Base Defense + Hero's Defense Skill) * Defend Stance Multiplier
    total_defense = (defender_base_defense + hero_defense_skill) * defend_stance_multiplier

    # Step 5: Calculate the damage bonus from the difference between Attack and Defense
    # Since Attack > Defense, the bonus is 5% for each point of difference.
    # Damage Bonus = (Total Attack - Total Defense) * 0.05
    attack_vs_defense_bonus = (total_attack - total_defense) * 0.05
    attack_vs_defense_multiplier = 1.0 + attack_vs_defense_bonus

    # Step 6: Calculate the final damage by applying all modifiers
    base_stack_damage = num_attackers * attacker_base_damage
    final_damage_unfloored = (
        base_stack_damage
        * attack_vs_defense_multiplier
        * racial_bonus_multiplier
        * spell_reduction_multiplier
    )
    final_damage = math.floor(final_damage_unfloored)
    
    # Step 7: Print the detailed calculation
    print("--- Damage Calculation Breakdown ---")
    print(f"Attacker's Total Attack: {attacker_base_attack} (base) + {hero_attack_skill} (hero) = {total_attack}")
    print(f"Defender's Total Defense: ({defender_base_defense} (base) + {hero_defense_skill} (hero)) * {defend_stance_multiplier} (defend) = {total_defense:.1f}")
    print(f"Attack vs. Defense Multiplier: 1 + ({total_attack} - {total_defense:.1f}) * 0.05 = {attack_vs_defense_multiplier:.2f}")

    print("\n--- Final Damage Equation ---")
    print(f"Total Damage = floor( [Base Stack Damage] * [Atk/Def Multiplier] * [Racial Bonus] * [Spell Reduction] )")
    print(f"Total Damage = floor( ({num_attackers} * {attacker_base_damage}) * {attack_vs_defense_multiplier:.2f} * {racial_bonus_multiplier:.1f} * {spell_reduction_multiplier:.1f} )")
    print(f"Total Damage = floor( {base_stack_damage} * {attack_vs_defense_multiplier:.2f} * {racial_bonus_multiplier:.1f} * {spell_reduction_multiplier:.1f} )")
    print(f"Total Damage = floor({final_damage_unfloored:.3f})")
    print(f"\nFinal calculated damage: {final_damage}")

calculate_homm3_damage()