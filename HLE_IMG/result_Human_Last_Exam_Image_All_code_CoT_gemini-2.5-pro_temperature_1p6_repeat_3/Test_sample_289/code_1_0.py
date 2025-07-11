# The elevation of Descanso, CA according to the Geographic Names Information System (GNIS)
elevation_feet = 3524

# The conversion factor from feet to meters
feet_to_meters_conversion = 0.3048

# Calculate the elevation in meters
elevation_meters = elevation_feet * feet_to_meters_conversion

# Print the full calculation and the result
print(f"The elevation of Descanso is {elevation_feet} feet.")
print(f"The conversion is: {elevation_feet} feet * {feet_to_meters_conversion} meters/foot = {elevation_meters:.2f} meters.")