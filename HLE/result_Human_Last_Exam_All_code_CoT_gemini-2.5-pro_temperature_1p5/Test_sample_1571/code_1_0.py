# This script calculates the total damage based on the optimal play described.

# --- Initial Damage ---
# Player A uses the Channel ability of Twinshot Sniper.
channel_damage = 2
print(f"Player A deals {channel_damage} damage with Twinshot Sniper's Channel ability.")

# --- Combat Damage Calculation ---
# Player A attacks with a board strengthened by pre-combat plays.
# Player B makes the optimal block to minimize damage taken.
# The following creatures are unblocked and deal combat damage:
specialist_damage = 4
boar_damage = 4
welder_damage = 3
apprentice_damage = 2

# Calculate the total combat damage
combat_damage = specialist_damage + boar_damage + welder_damage + apprentice_damage
print("The unblocked creatures deal the following damage:")
print(f"Replication Specialist: {specialist_damage}")
print(f"Ironhoof Boar: {boar_damage}")
print(f"Scrap Welder: {welder_damage}")
print(f"Iron Apprentice: {apprentice_damage}")

print(f"\nThe total combat damage is calculated by the equation: {specialist_damage} + {boar_damage} + {welder_damage} + {apprentice_damage} = {combat_damage}")

# --- Total Life Loss Calculation ---
total_damage = channel_damage + combat_damage
print("\n--- Final Calculation ---")
print("The total life Player B loses is the sum of the initial Channel damage and the combat damage.")
print(f"Final Equation: {channel_damage} + {combat_damage} = {total_damage}")