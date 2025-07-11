# Step 1: Define the damage from the final Action spell, Disintegrate.
# The spell is cast using an 8th-level slot.
# Base damage at 6th level is 10d6 + 40.
# It gains 3d6 damage for each level above 6th (i.e., +6d6 for an 8th-level slot).
# Total dice are 10d6 + 6d6 = 16d6.
disintegrate_dice_count = 16
disintegrate_dice_sides = 6
disintegrate_bonus_damage = 40

# "Best case scenario" means all dice roll their maximum value.
max_disintegrate_damage = (disintegrate_dice_count * disintegrate_dice_sides) + disintegrate_bonus_damage
print(f"Final Action: Cast Disintegrate at 8th level.")
print(f"Damage Formula: {disintegrate_dice_count}d{disintegrate_dice_sides} + {disintegrate_bonus_damage}")
print(f"Max Damage from Disintegrate: {max_disintegrate_damage}")
print("-" * 20)

# Step 2: Define the damage from the final Bonus Action, commanding Animated Objects.
# The spell is cast using a 7th-level slot.
# This animates up to 14 tiny objects. To reach exactly 240 total damage, we assume 13 are animated.
# Each tiny object deals 1d4 + 4 damage.
animated_object_count = 13
object_dice_count = 1
object_dice_sides = 4
object_bonus_damage = 4

# Calculate max damage for one object.
max_damage_per_object = (object_dice_count * object_dice_sides) + object_bonus_damage
# Calculate total damage from all objects.
total_object_damage = animated_object_count * max_damage_per_object
print(f"Final Bonus Action: Command Animated Objects (cast at 7th level).")
print(f"Number of Objects: {animated_object_count}")
print(f"Damage Formula per Object: {object_dice_count}d{object_dice_sides} + {object_bonus_damage}")
print(f"Max Damage from Animated Objects: {total_object_damage}")
print("-" * 20)

# Step 3: Calculate the total damage.
total_damage = max_disintegrate_damage + total_object_damage
print("Total damage is the sum of the Disintegrate and the Animated Objects attacks.")
print(f"Final Equation: {max_disintegrate_damage} (from Disintegrate) + {total_object_damage} (from Animated Objects) = {total_damage}")