# Step 1 & 2: Calculate pre-combat damage from Twinshot Sniper ETB triggers.
twinshot_sniper_etb_damage = 2
replicated_twinshot_sniper_etb_damage = 2
pre_combat_damage = twinshot_sniper_etb_damage + replicated_twinshot_sniper_etb_damage

# Step 3: Calculate combat damage.
# Player B must use their two blockers on the 3/4 Mukotai Soulripper due to Menace.
# The other attackers are unblocked.
ironhoof_boar_power = 4
replication_specialist_power = 2
scrap_welder_power = 1
twinshot_sniper_power = 2
twinshot_sniper_token_power = 2

combat_damage = (ironhoof_boar_power + 
                 replication_specialist_power + 
                 scrap_welder_power + 
                 twinshot_sniper_power + 
                 twinshot_sniper_token_power)

# Step 4: Calculate total life loss for Player B.
total_life_loss = pre_combat_damage + combat_damage

# Final Output: Print the equation showing how the total life loss is calculated.
print(f"{twinshot_sniper_etb_damage} + {replicated_twinshot_sniper_etb_damage} + {ironhoof_boar_power} + {replication_specialist_power} + {scrap_welder_power} + {twinshot_sniper_power} + {twinshot_sniper_token_power} = {total_life_loss}")