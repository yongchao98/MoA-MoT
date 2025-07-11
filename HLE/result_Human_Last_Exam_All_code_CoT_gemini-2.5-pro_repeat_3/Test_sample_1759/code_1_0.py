import math

# Step 1: Define all the base stats and conditions
num_archangels = 3
archangel_damage = 50
archangel_attack = 30
devil_defense = 21

red_hero_attack = 19
blue_hero_defense = 1

# Modifiers
defensive_stance_mod = 1.2  # +20%
racial_bonus_mod = 0.5  # +50%
spell_reduction_mod = 0.3  # -30%

# Step 2: Calculate effective attack and defense stats
total_attack = archangel_attack + red_hero_attack
# The defensive stance bonus creates a fractional defense value
effective_defense = (devil_defense + blue_hero_defense) * defensive_stance_mod

# Step 3: Calculate base damage for the stack
base_damage = num_archangels * archangel_damage

# Step 4: Calculate damage after Attack/Defense modifier
# Damage bonus = (Total Attack - Effective Defense) * 5%
ad_multiplier = 1 + (total_attack - effective_defense) * 0.05
damage_after_ad = math.floor(base_damage * ad_multiplier)

# Step 5: Add the racial bonus (+50% vs Devils)
# The bonus is calculated on the current damage, floored, and then added.
racial_bonus_amount = math.floor(damage_after_ad * racial_bonus_mod)
damage_after_racial = damage_after_ad + racial_bonus_amount

# Step 6: Subtract the spell reduction (Protection from Water, Expert)
# The reduction is calculated on the current damage, floored, and then subtracted.
spell_reduction_amount = math.floor(damage_after_racial * spell_reduction_mod)
final_damage = damage_after_racial - spell_reduction_amount

# Step 7: Print the detailed calculation steps
print("--- Damage Calculation Breakdown ---")
print(f"Base Damage: {num_archangels} Archangels * {archangel_damage} Damage = {base_damage}")
print(f"Total Attack: {archangel_attack} (Unit) + {red_hero_attack} (Hero) = {total_attack}")
print(f"Effective Defense: ({devil_defense} (Unit) + {blue_hero_defense} (Hero)) * {defensive_stance_mod} (Defend) = {effective_defense}")
print("\n--- Step-by-Step Damage ---")
print(f"1. Damage after A/D modifier: floor({base_damage} * (1 + ({total_attack} - {effective_defense}) * 0.05)) = {damage_after_ad}")
print(f"2. Damage after Racial Bonus: {damage_after_ad} + floor({damage_after_ad} * {racial_bonus_mod}) = {damage_after_ad} + {racial_bonus_amount} = {damage_after_racial}")
print(f"3. Final Damage after Spell Reduction: {damage_after_racial} - floor({damage_after_racial} * {spell_reduction_mod}) = {damage_after_racial} - {spell_reduction_amount} = {final_damage}")

print(f"\nFinal calculated damage: {final_damage}")
<<<A>>>