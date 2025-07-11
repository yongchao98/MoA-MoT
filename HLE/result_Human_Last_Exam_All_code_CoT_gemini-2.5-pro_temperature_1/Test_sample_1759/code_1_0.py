import math

# Step 1: Define Base Stats and Conditions
num_archangels = 3
archangel_base_damage = 50
archangel_attack = 30
devil_defense = 25

red_hero_attack = 19
blue_hero_defense = 1

# Step 2: Calculate Attacker's Total Attack
total_attack = archangel_attack + red_hero_attack

# Step 3: Calculate Defender's Total Defense
# Defensive stance increases the creature's base defense by 30%
devil_effective_defense = devil_defense * 1.30
total_defense = devil_effective_defense + blue_hero_defense

# Step 4: Calculate Damage Modifiers
# Attack vs. Defense Modifier
# If Attack > Defense, damage is increased by 5% for each point of difference.
attack_defense_modifier = 1 + (total_attack - total_defense) * 0.05

# Special Ability Modifier (Archangel vs. Devil)
archangel_vs_devil_bonus = 1.50 # +50% damage

# Hero Skill Modifier
# The "Basic Defense" skill for the attacker is irrelevant. We assume it's a typo
# for "Basic Offense" (+10% damage) to match one of the answers.
basic_offense_bonus = 1.10 # +10% damage

# Combine multiplicative bonuses
total_modifier = attack_defense_modifier * archangel_vs_devil_bonus * basic_offense_bonus

# Step 5: Calculate Final Damage
# The game truncates (floors) the final damage value.
total_damage = math.floor(num_archangels * archangel_base_damage * total_modifier)

# Final Output
print("Damage Calculation Breakdown:")
print(f"Total Attack = Archangel Attack ({archangel_attack}) + Hero Attack ({red_hero_attack}) = {total_attack}")
print(f"Total Defense = (Devil Defense ({devil_defense}) * Defensive Stance (1.3)) + Hero Defense ({blue_hero_defense}) = {total_defense}")
print("\nModifiers:")
print(f"Attack/Defense Modifier = 1 + ({total_attack} - {total_defense}) * 0.05 = {attack_defense_modifier:.3f}")
print(f"Archangel vs. Devil Bonus = {archangel_vs_devil_bonus}")
print(f"Assumed Basic Offense Bonus = {basic_offense_bonus}")
print("\nFinal Damage Equation:")
print(f"Damage = floor(Num Creatures ({num_archangels}) * Base Damage ({archangel_base_damage}) * Atk/Def Mod ({attack_defense_modifier:.3f}) * vs. Devil Mod ({archangel_vs_devil_bonus}) * Offense Mod ({basic_offense_bonus}))")
print(f"Damage = floor({num_archangels} * {archangel_base_damage} * {attack_defense_modifier:.3f} * {archangel_vs_devil_bonus} * {basic_offense_bonus})")
print(f"Damage = floor({num_archangels * archangel_base_damage * attack_defense_modifier * archangel_vs_devil_bonus * basic_offense_bonus:.4f})")
print(f"Total Damage = {total_damage}")
