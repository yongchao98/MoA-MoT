# Given data
laps = 7637264
distance_per_lap_km = 0.25
water_per_km_ml = 60

# Calculate total distance
total_distance_km = laps * distance_per_lap_km

# Calculate total water needed
total_water_ml = total_distance_km * water_per_km_ml

# Output the result
print(total_water_ml)