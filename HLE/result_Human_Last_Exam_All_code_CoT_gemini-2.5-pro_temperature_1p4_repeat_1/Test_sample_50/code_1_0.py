# Spell slot levels used in this combo
slot_level_dbf = 7
slot_level_ofs = 6
slot_level_ahw = 8

# Delayed Blast Fireball (7th-level)
dbf_base_dice = 12
# The spell was cast on turn 1 and held until turn 3, so it gains power over 2 full rounds.
rounds_held = 2
dbf_bonus_dice = rounds_held
dbf_total_dice = dbf_base_dice + dbf_bonus_dice
# A d6 die is used for damage
dbf_die_type = 6
dbf_max_damage = dbf_total_dice * dbf_die_type
print(f"Delayed Blast Fireball ({dbf_total_dice}d{dbf_die_type}) damage: {dbf_max_damage}")

# Otiluke's Freezing Sphere (6th-level)
ofs_dice = 10
# A d6 die is used for damage
ofs_die_type = 6
ofs_max_damage = ofs_dice * ofs_die_type
print(f"Otiluke's Freezing Sphere ({ofs_dice}d{ofs_die_type}) damage: {ofs_max_damage}")

# Abi-Dalzim's Horrid Wilting (8th-level)
ahw_dice = 12
# A d8 die is used for damage
ahw_die_type = 8
ahw_max_damage = ahw_dice * ahw_die_type
print(f"Abi-Dalzim's Horrid Wilting ({ahw_dice}d{ahw_die_type}) damage: {ahw_max_damage}")

# Total damage calculation
total_damage = dbf_max_damage + ofs_max_damage + ahw_max_damage
print(f"\nFinal Equation: {dbf_max_damage} (Delayed Blast Fireball) + {ofs_max_damage} (Otiluke's Freezing Sphere) + {ahw_max_damage} (Abi-Dalzim's Horrid Wilting)")
print(f"Total Maximum Damage: {total_damage}")