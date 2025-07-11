import math

# --- 1. Define Base Stats and Conditions ---
# Attacker: Red Player (Thant)
num_archangels = 3
archangel_damage = 50
archangel_attack = 30
red_hero_attack = 19

# Defender: Blue Player (Bron)
devil_defense = 21
blue_hero_defense = 1
defensive_stance_bonus = 1.2 # +20% defense

# Special bonus
archangel_vs_devil_bonus = 1.5 # +50% damage

# --- 2. Calculate Total Attack Skill ---
total_attack_skill = archangel_attack + red_hero_attack

# --- 3. Calculate Total Defense Skill ---
# The defensive stance bonus is applied to the combined creature and hero defense.
base_defense_for_calc = devil_defense + blue_hero_defense
total_defense_skill = base_defense_for_calc * defensive_stance_bonus

# --- 4. Calculate Base Damage ---
base_damage = num_archangels * archangel_damage

# --- 5. Apply Damage Formula with Skill Difference ---
skill_difference = total_attack_skill - total_defense_skill

# Damage is increased by 5% for each point of attack over defense
attack_over_defense_multiplier = 1 + (skill_difference * 0.05)
damage_after_skills = math.floor(base_damage * attack_over_defense_multiplier)

# --- 6. Apply Racial Bonus ---
final_damage = math.floor(damage_after_skills * archangel_vs_devil_bonus)

# --- 7. Print the step-by-step calculation ---
print("Calculating Damage: 3 Archangels vs. 3 Devils")
print("-" * 40)
print(f"Total Attack Skill = {archangel_attack} (Archangel) + {red_hero_attack} (Hero) = {total_attack_skill}")
print(f"Total Defense Skill = ({devil_defense} (Devil) + {blue_hero_defense} (Hero)) * {defensive_stance_bonus} (Defend Stance) = {total_defense_skill:.1f}")
print("-" * 40)
print(f"Base Damage for Stack = {num_archangels} Archangels * {archangel_damage} damage = {base_damage}")
print(f"Damage after A/D bonus = floor({base_damage} * (1 + ({total_attack_skill} - {total_defense_skill:.1f}) * 0.05)) = {damage_after_skills}")
print(f"Final Damage after Racial Bonus = floor({damage_after_skills} * {archangel_vs_devil_bonus}) = {final_damage}")
print("-" * 40)
print(f"The final calculated damage is {final_damage}.")
