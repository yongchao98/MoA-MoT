# Step 1: Define the power of each creature that will attack.
# Player A clears the blockers and creates two Iron Apprentice creatures before combat.
replication_specialist_power = 2
mukotai_soulripper_power = 4
scrap_welder_power = 3
ironhoof_boar_power = 4
iron_apprentice_1_power = 1
iron_apprentice_2_power = 1

# Step 2: Calculate the total combat damage.
combat_damage = (replication_specialist_power + 
                 mukotai_soulripper_power + 
                 scrap_welder_power + 
                 ironhoof_boar_power + 
                 iron_apprentice_1_power + 
                 iron_apprentice_2_power)

# Step 3: Define the additional damage from Twinshot Sniper's Channel ability.
sniper_channel_damage = 2

# Step 4: Calculate the total damage dealt to Player B.
total_damage = combat_damage + sniper_channel_damage

# Step 5: Print the final equation and the result.
print(f"{replication_specialist_power} (Replication Specialist) + "
      f"{mukotai_soulripper_power} (Mukotai Soulripper) + "
      f"{scrap_welder_power} (Scrap Welder) + "
      f"{ironhoof_boar_power} (Ironhoof Boar) + "
      f"{iron_apprentice_1_power} (Iron Apprentice) + "
      f"{iron_apprentice_2_power} (Iron Apprentice Token) + "
      f"{sniper_channel_damage} (Twinshot Sniper) = {total_damage}")
