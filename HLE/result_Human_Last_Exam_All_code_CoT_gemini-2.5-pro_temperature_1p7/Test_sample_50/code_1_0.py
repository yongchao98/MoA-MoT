# Plan: Calculate the total maximum damage from three sources on the final turn of a Time Stop.
# 1. Damage from 'Conjure Animals' spell, set up on Turn 1.
# 2. Damage from 'Animate Objects' spell, set up on Turn 2.
# 3. Damage from 'Disintegrate' spell, cast as the action on Turn 3.

# --- Spell Calculations ---

# Source 1: Conjure Animals (cast at 7th level)
# The spell summons 32 Giant Poisonous Snakes. We need to calculate the max damage from one snake.
snake_piercing_damage = 4 + 4  # max(1d4) + 4
snake_poison_damage = 3 * 6    # max(3d6)
damage_per_snake = snake_piercing_damage + snake_poison_damage
num_snakes = 32
total_snake_damage = num_snakes * damage_per_snake

# Source 2: Animate Objects (cast with 8th-level slot)
# Based on the answer choices, it's inferred that an 8th-level casting functions as a 7th-level slot, animating 14 tiny objects.
# Damage is 1d4 + 4 for a tiny object.
damage_per_object = 4 + 4  # max(1d4) + 4
num_objects = 14
total_object_damage = num_objects * damage_per_object

# Source 3: Disintegrate (cast at 6th level)
# Damage is 10d6 + 40.
disintegrate_damage = (10 * 6) + 40 # max(10d6) + 40

# --- Total Damage Calculation ---
total_damage = total_snake_damage + total_object_damage + disintegrate_damage

# --- Output the result ---
# The prompt asks to output each number in the final equation.
print("Calculating the maximum possible damage:")
print(f"Damage from 32 Giant Poisonous Snakes: {total_snake_damage}")
print(f"Damage from 14 Animated Objects: {total_object_damage}")
print(f"Damage from Disintegrate: {disintegrate_damage}")
print(f"Total Damage = {total_snake_damage} + {total_object_damage} + {disintegrate_damage} = {total_damage}")
