import math

# Step 1: Define base stats and conditions
num_attackers = 3
attacker_base_damage = 50  # Archangel fixed damage
attacker_base_attack = 30  # Archangel base attack
defender_base_defense = 21 # Devil base defense

attacker_hero_offense = 19 # Red Hero's Offense
defender_hero_defense = 1  # Blue Hero's Defense

# Step 2: Calculate effective attack and defense stats
effective_attack = attacker_base_attack + attacker_hero_offense

# Defender's defense is increased by "Defensive Stance" (20%)
# This bonus applies to the combined unit and hero defense
defender_pre_stance_defense = defender_base_defense + defender_hero_defense
defensive_stance_multiplier = 1.2
effective_defense = defender_pre_stance_defense * defensive_stance_multiplier

# Step 3: Calculate the damage bonus from the attack/defense difference
# Bonus is 5% for each point of Attack over Defense
attack_advantage = effective_attack - effective_defense
damage_multiplier = 1 + (attack_advantage * 0.05)

# Step 4: Calculate damage before racial bonus, then truncate
damage_before_racial_bonus = num_attackers * attacker_base_damage * damage_multiplier
truncated_damage = math.floor(damage_before_racial_bonus)

# Step 5: Apply the racial bonus (+50% vs Devils) and truncate the final result
racial_bonus_multiplier = 1.5
final_damage = truncated_damage * racial_bonus_multiplier
final_damage_truncated = math.floor(final_damage)

# --- Output the breakdown ---
print("--- Damage Calculation Breakdown ---")
print(f"1. Attacker's Effective Attack: {attacker_base_attack} (Archangel) + {attacker_hero_offense} (Hero) = {effective_attack}")
print(f"2. Defender's Defense before Stance: {defender_base_defense} (Devil) + {defender_hero_defense} (Hero) = {defender_pre_stance_defense}")
print(f"3. Defender's Effective Defense with Stance: {defender_pre_stance_defense} * {defensive_stance_multiplier} (20% bonus) = {effective_defense:.1f}")
print(f"4. Skill-Based Damage Multiplier: 1 + ({effective_attack} - {effective_defense:.1f}) * 0.05 = {damage_multiplier:.2f}")
print(f"5. Damage before Racial Bonus: {num_attackers} * {attacker_base_damage} * {damage_multiplier:.2f} = {damage_before_racial_bonus:.2f}")
print(f"6. Truncated Damage: floor({damage_before_racial_bonus:.2f}) = {truncated_damage}")
print(f"7. Applying +50% Racial Bonus vs. Devils: {truncated_damage} * {racial_bonus_multiplier} = {final_damage:.1f}")
print(f"8. Final Damage (Truncated): floor({final_damage:.1f}) = {final_damage_truncated}")
print("\n--- Final Answer ---")
print(f"The total damage inflicted by 3 Archangels on 3 Devils is {final_damage_truncated}.")

# --- Output the final equation ---
print("\nFinal Equation:")
# Using variables for clarity in the equation string
A_base = attacker_base_attack
H_off = attacker_hero_offense
D_base = defender_base_defense
H_def = defender_hero_defense
N_att = num_attackers
D_dmg = attacker_base_damage
stance_mult = defensive_stance_multiplier
hate_mult = racial_bonus_multiplier
final_val = final_damage_truncated

print(f"floor( floor( {N_att} * {D_dmg} * (1 + ({A_base} + {H_off} - ({D_base} + {H_def}) * {stance_mult}) * 0.05) ) * {hate_mult} ) = {final_val}")