# Player A's initial attacking creatures' power
replication_specialist_power = 3
ironhoof_boar_power = 4
ninja_token_1_power = 1
ninja_token_2_power = 1

# Step 1: Use Scrap Welder and Iron Apprentice
# No change to the final list of attackers, but this creates a 4/4 Iron Apprentice.
# This 4/4 creature will be used to crew the vehicle.

# Step 2: Remove blockers
# Twinshot Sniper's Channel ability removes one 1/1 blocker.
# March of Wretched Sorrow removes the other 1/1 blocker.
# Player B now has no blockers.

# Step 3: Crew the vehicle
# The 4/4 Iron Apprentice is tapped to crew the Mukotai Soulripper.
# Mukotai Soulripper becomes a 4/4 creature.
mukotai_soulripper_initial_power = 4

# Step 4: Declare attackers and use Mukotai Soulripper's ability
# The attacking creatures are:
# - Mukotai Soulripper (4/4)
# - Replication Specialist (3/4)
# - Ironhoof Boar (4/3)
# - Two Ninja tokens (1/1 each)
#
# When Mukotai Soulripper attacks, its ability triggers.
# We sacrifice the tapped 4/4 Iron Apprentice (which was used for crew and is not attacking).
# This puts a +1/+1 counter on Mukotai Soulripper.
mukotai_soulripper_final_power = mukotai_soulripper_initial_power + 1

# Step 5: Calculate the total damage dealt by the unblocked attackers.
# The attackers are Mukotai Soulripper, Replication Specialist, Ironhoof Boar, and two Ninjas.
total_life_loss = (
    mukotai_soulripper_final_power +
    replication_specialist_power +
    ironhoof_boar_power +
    ninja_token_1_power +
    ninja_token_2_power
)

# Output the final calculation
print("Player B's total life loss is calculated as follows:")
print(f"Damage from Mukotai Soulripper: {mukotai_soulripper_final_power}")
print(f"Damage from Replication Specialist: {replication_specialist_power}")
print(f"Damage from Ironhoof Boar: {ironhoof_boar_power}")
print(f"Damage from first Ninja token: {ninja_token_1_power}")
print(f"Damage from second Ninja token: {ninja_token_2_power}")
print(f"Final Equation: {mukotai_soulripper_final_power} + {replication_specialist_power} + {ironhoof_boar_power} + {ninja_token_1_power} + {ninja_token_2_power}")
print(f"Total life loss for Player B: {total_life_loss}")
