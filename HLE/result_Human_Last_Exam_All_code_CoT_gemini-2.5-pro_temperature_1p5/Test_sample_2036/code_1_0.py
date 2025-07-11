# ATC Separation Calculation for Bilbao Approach

# Step 1: Define the constants based on ICAO standards and common ATC procedures.
# The problem involves a circling approach. The protected area for this maneuver is key.
# We assume a Category C aircraft (like an A320/B737), which is standard for this type of airport.
# ICAO PANS-OPS (Doc 8168) specifies a circling radius of 4.2 NM for a Category C aircraft.
icao_circling_radius_cat_c = 4.2  # Nautical Miles

# Standard radar separation used in approach control is 5 NM.
# This is the buffer needed between the departing aircraft and the circling aircraft's protected area.
standard_radar_separation = 5.0  # Nautical Miles

# Step 2: Calculate the total required distance from the VOR.
# The VOR is on the airfield, serving as the center of the conflict area.
# To clear a departure, the arriving aircraft must be far enough away to allow the
# departure to take off and climb without infringing on the 5 NM separation standard
# from the arriving aircraft's protected area.
# Therefore, the distance is the sum of the circling radius and the separation standard.
total_distance_from_vor = icao_circling_radius_cat_c + standard_radar_separation

# Step 3: Print the result and the equation used for clarity.
# The user needs to see the final answer and how it was derived.
print(f"To provide safe separation, the arriving traffic must be at least {total_distance_from_vor:.1f} NM from the VOR.")
print("This distance is the sum of the aircraft's protected circling area radius and the standard radar separation.")
print("The final equation is:")
print(f"{total_distance_from_vor:.1f} NM = {icao_circling_radius_cat_c} NM (Circling Radius) + {standard_radar_separation} NM (Separation Standard)")