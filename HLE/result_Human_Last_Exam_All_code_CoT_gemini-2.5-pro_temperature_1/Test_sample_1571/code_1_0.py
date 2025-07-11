# Player A's goal is to maximize the damage dealt to Player B.
# This is achieved by spending mana on creatures and abilities rather than direct damage from hand.

# Step 1: Calculate the power of each creature after Player A's main phase actions.
# Ironhoof Boar is enhanced by Mukotai Soulripper.
# Base power of Ironhoof Boar is 4.
# It gets +2 power from Mukotai Soulripper being equipped.
# It gets another +2 power from Mukotai's triggered ability when a 1/1 Ninja is sacrificed.
boar_damage = 4 + 2 + 2

# The power of other attacking creatures.
specialist_damage = 2
welder_damage = 3
# Player A creates two 2/2 Iron Apprentice creatures.
apprentice_damage = 2 + 2
# One 1/1 Ninja token attacks, the other was sacrificed.
ninja_damage = 1

# Step 2: Sum the damage from all unblocked attackers.
# Player B has only flying blockers, which are ineffective against Player A's ground creatures.
total_damage_dealt = boar_damage + specialist_damage + welder_damage + apprentice_damage + ninja_damage

# Step 3: Print the breakdown of the damage calculation.
print("To maximize damage, Player A takes the following steps:")
print("1. Casts Iron Apprentice, creating a copy with Replication Specialist.")
print("2. Reconfigures Mukotai Soulripper onto Ironhoof Boar.")
print("3. In combat, sacrifices a Ninja token to boost the Boar's power.")
print("\nPlayer B's flying blockers cannot stop the ground assault. The total life lost is calculated as follows:")
print(f"Damage from Ironhoof Boar: {boar_damage}")
print(f"Damage from Replication Specialist: {specialist_damage}")
print(f"Damage from Scrap Welder: {welder_damage}")
print(f"Damage from two Iron Apprentices: {apprentice_damage}")
print(f"Damage from one Ninja token: {ninja_damage}")
print(f"\nTotal life lost by Player B: {boar_damage} + {specialist_damage} + {welder_damage} + {apprentice_damage} + {ninja_damage} = {total_damage_dealt}")