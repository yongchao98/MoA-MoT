# Step 1: Define the damage sources and spell details based on the plan.

# Damage from Giant Poisonous Snake (CR 1/4)
# Bite: 1d4 piercing + 4
# Poison: 3d6 poison on a failed DC 11 Con save.
# Best case: Max roll on dice, target auto-fails save.
snake_piercing_damage = 4 + 4  # max roll on 1d4 is 4
snake_poison_damage = 3 * 6   # max roll on 3d6
damage_per_snake = snake_piercing_damage + snake_poison_damage

# Step 2: Calculate damage from summoned creatures.
# 'Conjure Animals' is cast with an 8th-level slot.
# Assuming generous upcasting, this works like a 9th-level slot.
# A 9th-level slot summons 4x the base number of CR 1/4 beasts.
# Base number is 8.
num_snakes = 4 * 8
total_snake_damage = num_snakes * damage_per_snake

# Step 3: Calculate damage from the final turn's spells.
# Remaining slots: 7, 6, 5, 4, 3, 2, 1.
# Action: Scorching Ray with the 7th-level slot.
# Creates 3 base rays + (7-2) = 8 rays total.
# Each ray is 2d6 fire. Crit damage is 4d6.
sr_rays = 3 + (7 - 2)
sr_damage_per_ray = 4 * 6 # max roll on 4d6 (crit)
total_sr_damage = sr_rays * sr_damage_per_ray

# Bonus Action: Spiritual Weapon with the 6th-level slot.
# Base damage 1d8. Upcast: +1d8 for every 2 levels above 2nd (4th, 6th).
# 6th-level is +2d8. Total damage is 3d8 force. Crit damage is 6d8.
sw_damage_dice = 6 * 8 # max roll on 6d8 (crit)
total_sw_damage = sw_damage_dice

# Step 4: Sum all damage sources and print the result.
total_damage = total_snake_damage + total_sr_damage + total_sw_damage

print(f"The plan is to summon {num_snakes} Giant Poisonous Snakes, then attack with Scorching Ray and Spiritual Weapon.")
print(f"Maximized damage from {num_snakes} snakes: {num_snakes} * ({snake_piercing_damage} + {snake_poison_damage}) = {total_snake_damage}")
print(f"Maximized damage from a critical 7th-level Scorching Ray ({sr_rays} rays): {sr_rays} * {sr_damage_per_ray} = {total_sr_damage}")
print(f"Maximized damage from a critical 6th-level Spiritual Weapon: {total_sw_damage}")
print("\nFinal Equation:")
print(f"{total_snake_damage} (Summons) + {total_sr_damage} (Scorching Ray) + {total_sw_damage} (Spiritual Weapon) = {total_damage}")
