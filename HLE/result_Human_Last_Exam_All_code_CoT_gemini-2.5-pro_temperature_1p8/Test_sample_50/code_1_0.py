# Step 1: Define parameters and spell damage based on the plan.
# Assume best-case scenario for all dice rolls (maximum value).
# Character has +5 spellcasting modifier, but it doesn't apply to these damage rolls.

# Spell 1: Crown of Stars (cast on Turn 1)
# 7th-level spell. No concentration.
# Creates 7 motes. Each deals 4d12 damage.
motes_count = 7
damage_die_crown = 12
dice_per_mote = 4
max_damage_per_mote = dice_per_mote * damage_die_crown
total_max_damage_crown = motes_count * max_damage_per_mote

print("Turn 1: Cast Crown of Stars (7th-level slot).")
print("Turn 2: Cast Bigby's Hand (using 8th-level slot).")
print("\n--- Final Turn ---")

# On Turn 3, use Bonus Action to fire all motes.
print(f"Bonus Action: Crown of Stars attacks with {motes_count} motes.")
print(f"Damage per mote is {dice_per_mote}d{damage_die_crown}. Max damage: {dice_per_mote} * {damage_die_crown} = {max_damage_per_mote}")
print(f"Total Crown of Stars Damage: {motes_count} * {max_damage_per_mote} = {total_max_damage_crown}\n")

# Spell 2: Bigby's Hand (cast on Turn 2, used on Turn 3)
# 5th-level spell, upcast using an 8th-level slot.
# Base damage (Clenched Fist) is 4d8 force damage.
# Upcasting adds 2d8 damage for each slot level above 5th.
# 8th level is 3 levels above 5th, so add 3 * 2d8 = 6d8.
# Total damage is 4d8 + 6d8 = 10d8.
damage_die_hand = 8
base_dice_hand = 4
upcast_level = 8
base_level_hand = 5
extra_dice_per_level = 2
added_dice_hand = (upcast_level - base_level_hand) * extra_dice_per_level
total_dice_hand = base_dice_hand + added_dice_hand
max_damage_hand = total_dice_hand * damage_die_hand

# On Turn 3, use Action for Clenched Fist attack.
print(f"Action: Bigby's Hand (upcast to level {upcast_level}) attacks.")
print(f"Damage is {total_dice_hand}d{damage_die_hand}. Max damage: {total_dice_hand} * {damage_die_hand} = {max_damage_hand}\n")

# Step 2: Calculate total maximum damage.
total_damage = total_max_damage_crown + max_damage_hand

print("--- Total Damage Calculation ---")
print(f"Final Equation: ({motes_count} * {dice_per_mote} * {damage_die_crown}) + ({total_dice_hand} * {damage_die_hand}) = {total_max_damage_crown} + {max_damage_hand}")
print(f"Total Maximum Damage: {total_damage}")