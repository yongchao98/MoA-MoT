# This script calculates the maximum life Player B will lose in a single turn
# based on the provided Magic: The Gathering scenario.

# --- Damage Calculation ---

# 1. Combat Damage from Ironhoof Boar
# Ironhoof Boar has a base power of 4 and gets +2/+0 when it attacks.
boar_base_power = 4
boar_attack_bonus = 2
boar_damage = boar_base_power + boar_attack_bonus

# 2. Combat Damage from Mukotai Soulripper
# Mukotai Soulripper has a base power of 4.
# Player A sacrifices the Iron Apprentice token to its attack trigger. This results in two bonuses:
# - The Soulripper's own ability gives it +2/+0.
# - The sacrificed Apprentice's death trigger puts its one +1/+1 counter on the Soulripper.
soulripper_base_power = 4
soulripper_attack_bonus = 2
apprentice_counter_bonus = 1
soulripper_damage = soulripper_base_power + soulripper_attack_bonus + apprentice_counter_bonus

# 3. Combat Damage from Replication Specialist
# Replication Specialist attacks for its base power.
specialist_damage = 2

# 4. Combat Damage from Iron Apprentice
# The original Iron Apprentice attacks for its power.
apprentice_damage = 1

# 5. Post-Combat Damage from Twinshot Sniper
# Player A uses the Channel ability of Twinshot Sniper.
sniper_channel_damage = 2

# --- Final Calculation and Output ---
total_life_lost = boar_damage + soulripper_damage + specialist_damage + apprentice_damage + sniper_channel_damage

print("The maximum life lost by Player B is calculated by summing damage from all sources:")
print(f"Ironhoof Boar Damage: {boar_damage}")
print(f"Mukotai Soulripper Damage: {soulripper_damage}")
print(f"Replication Specialist Damage: {specialist_damage}")
print(f"Iron Apprentice Damage: {apprentice_damage}")
print(f"Twinshot Sniper Channel Damage: {sniper_channel_damage}")
print("\nThe final equation is:")
print(f"{boar_damage} + {soulripper_damage} + {specialist_damage} + {apprentice_damage} + {sniper_channel_damage} = {total_life_lost}")
print(f"\nPlayer B will lose a total of {total_life_lost} life.")