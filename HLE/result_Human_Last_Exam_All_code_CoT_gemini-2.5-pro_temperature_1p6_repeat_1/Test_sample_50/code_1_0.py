# Character and spell parameters
spell_modifier = 5
num_turns = 3

# Turn 1: Cast Animate Objects (8th level)
# This animates 15 Tiny objects and 1 Small object.

# Turn 2: Hold action, maintain concentration.

# Turn 3: Final Burst
# Action: Scorching Ray (7th level)
# Bonus Action: Command Animated Objects

# Animate Objects Damage Calculation
# Best case roll = critical hit (dice are doubled) + max dice value
tiny_object_die_type = 4
tiny_object_num_dice = 1
num_tiny_objects = 15
tiny_object_crit_damage = (tiny_object_num_dice * 2 * tiny_object_die_type) + spell_modifier
total_tiny_object_damage = num_tiny_objects * tiny_object_crit_damage
print(f"1. Animate Objects (Tiny): {num_tiny_objects} objects deal {tiny_object_crit_damage} damage each (2d4+{spell_modifier} on a crit).")
print(f"   Max damage for one tiny object's critical hit: (4 + 4) + 5 = {tiny_object_crit_damage}")
print(f"   Equation: {num_tiny_objects} * {tiny_object_crit_damage} = {total_tiny_object_damage}")
print("-" * 20)

small_object_die_type = 8
small_object_num_dice = 1
num_small_objects = 1
small_object_crit_damage = (small_object_num_dice * 2 * small_object_die_type) + spell_modifier
total_small_object_damage = num_small_objects * small_object_crit_damage
print(f"2. Animate Objects (Small): {num_small_objects} object deals {small_object_crit_damage} damage (2d8+{spell_modifier} on a crit).")
print(f"   Max damage for one small object's critical hit: (8 + 8) + 5 = {small_object_crit_damage}")
print(f"   Equation: {num_small_objects} * {small_object_crit_damage} = {total_small_object_damage}")
print("-" * 20)

total_ao_damage = total_tiny_object_damage + total_small_object_damage
print(f"Total damage from all Animated Objects: {total_tiny_object_damage} + {total_small_object_damage} = {total_ao_damage}")
print("=" * 20)

# Scorching Ray Damage Calculation
# Cast at 7th level, it creates 8 rays. Each is a crit.
sr_level = 7
sr_base_rays = 3
sr_rays = sr_base_rays + (sr_level - 2)
sr_die_type = 6
sr_num_dice = 2
sr_crit_damage_per_ray = sr_num_dice * 2 * sr_die_type
total_sr_damage = sr_rays * sr_crit_damage_per_ray
print(f"3. Scorching Ray (7th level): Creates {sr_rays} rays. Each crits for {sr_crit_damage_per_ray} damage (4d6).")
print(f"   Max damage per ray's critical hit: 4 * 6 = {sr_crit_damage_per_ray}")
print(f"   Equation: {sr_rays} * {sr_crit_damage_per_ray} = {total_sr_damage}")
print("=" * 20)

# Final Total Damage Calculation
total_damage = total_ao_damage + total_sr_damage
print(f"Final Damage Calculation:")
print(f"Animate Objects ({total_ao_damage}) + Scorching Ray ({total_sr_damage}) = {total_damage}")