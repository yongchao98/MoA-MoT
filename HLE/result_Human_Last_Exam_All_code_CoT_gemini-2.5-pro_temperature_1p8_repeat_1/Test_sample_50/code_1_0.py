# Spell slot levels available: 1 through 8
# Time Stop allows for 3 turns.
# Damage will be dealt on the final turn using one Action and one Bonus Action.

# --- Bonus Action Damage: Animate Objects ---
# Spell is cast with the 8th-level slot.
# Casting at 5th level animates 10 objects.
# Upcasting adds 2 objects per level above 5th.
# Levels above 5th for an 8th-level slot: 8 - 5 = 3
# Additional objects: 3 * 2 = 6
# Total objects: 10 + 6 = 16

num_objects = 16
# Each tiny object deals 1d4 + 4 damage.
# Max damage on 1d4 is 4.
max_damage_per_object_die = 4
damage_modifier_per_object = 4
damage_per_object = max_damage_per_object_die + damage_modifier_per_object
total_damage_animate_objects = num_objects * damage_per_object

print(f"Damage from Animate Objects (8th level slot):")
print(f"{num_objects} objects * ({max_damage_per_object_die} + {damage_modifier_per_object}) damage per object = {total_damage_animate_objects}")

# --- Action Damage: Disintegrate ---
# The highest remaining spell slot is 7th level.
# Disintegrate (base 6th level) is cast with the 7th-level slot.
# Base damage: 10d6 + 40
# Upcast bonus for 7th level: +3d6
# Total dice: 13d6 + 40

num_dice_disintegrate = 13
max_roll_d6 = 6
damage_modifier_disintegrate = 40
damage_from_dice_disintegrate = num_dice_disintegrate * max_roll_d6
total_damage_disintegrate = damage_from_dice_disintegrate + damage_modifier_disintegrate

print(f"\nDamage from Disintegrate (7th level slot):")
print(f"({num_dice_disintegrate} * {max_roll_d6}) + {damage_modifier_disintegrate} = {total_damage_disintegrate}")


# --- Total Damage Calculation ---
total_damage = total_damage_animate_objects + total_damage_disintegrate

print(f"\nFinal Damage Calculation:")
print(f"Animate Objects Damage ({total_damage_animate_objects}) + Disintegrate Damage ({total_damage_disintegrate}) = {total_damage}")