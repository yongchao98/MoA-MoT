# Player A's plan to maximize damage to Player B.

# Initial state
player_b_life_loss = 0
scrap_welder_power = 3
mukotai_soulripper_power = 4
replication_specialist_power = 2
ironhoof_boar_power = 4

# --- Player A's First Main Phase ---

# 1. Player A uses the Channel ability of Twinshot Sniper.
# This costs {1}{R} and deals 2 damage directly to Player B.
twinshot_sniper_damage = 2
player_b_life_loss += twinshot_sniper_damage
print(f"Step 1: Player A channels Twinshot Sniper, dealing {twinshot_sniper_damage} direct damage.")
print(f"Player B has lost {player_b_life_loss} life so far.")
print("-" * 20)

# 2. Player A activates the ability of Clawing Torment on Scrap Welder.
# It costs {1}{B} to give +1/+1. Player A has {B}{B} and activates it twice.
activations = 2
scrap_welder_power += activations
print(f"Step 2: Player A activates Clawing Torment's ability {activations} times on Scrap Welder.")
print(f"Scrap Welder's power increases from 3 to {scrap_welder_power}.")
print("-" * 20)

# --- Player A's Combat Phase ---

# 3. Player A declares attackers. Mukotai is crewed by the two 1/1 Ninjas.
# Attackers: Replication Specialist (2), Scrap Welder (5), Ironhoof Boar (4), Mukotai Soulripper (4).
print("Step 3: Player A crews Mukotai Soulripper and declares attackers.")

# 4. Mukotai Soulripper's attack trigger resolves.
# Player A sacrifices a Ninja token to give Mukotai +1/+1.
mukotai_soulripper_power += 1
print(f"Step 4: Player A sacrifices a creature to Mukotai's trigger, making it a {mukotai_soulripper_power}/6 with Menace.")
print("-" * 20)

# 5. Player B declares blockers.
# With two 1/1 blockers, Player B must use both to block the Mukotai due to Menace.
print("Step 5: Player B must use their two 1/1 creatures to block the {mukotai_soulripper_power}/6 Mukotai with Menace.")
print("The other attackers are unblocked.")
print("-" * 20)

# 6. Combat damage is calculated.
# Damage from unblocked creatures: Replication Specialist + Scrap Welder + Ironhoof Boar
unblocked_combat_damage = replication_specialist_power + scrap_welder_power + ironhoof_boar_power
player_b_life_loss += unblocked_combat_damage
print("Step 6: Unblocked creatures deal combat damage.")
print(f"Combat Damage = {replication_specialist_power} (Replication Specialist) + {scrap_welder_power} (Scrap Welder) + {ironhoof_boar_power} (Ironhoof Boar) = {unblocked_combat_damage}")
print("-" * 20)

# --- Final Calculation ---
print("Final Calculation:")
print(f"Total Life Loss = Direct Damage ({twinshot_sniper_damage}) + Combat Damage ({unblocked_combat_damage})")
print(f"Total Life Loss = {twinshot_sniper_damage} + {unblocked_combat_damage} = {player_b_life_loss}")
