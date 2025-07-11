# Elevation of Descanso, CA in feet, according to GNIS
elevation_feet = 3524

# Conversion factor from feet to meters
feet_to_meters_conversion = 0.3048

# Calculate the elevation in meters
elevation_meters = elevation_feet * feet_to_meters_conversion

# Print the equation and the result
print(f"The elevation of Descanso is {elevation_feet} feet.")
print(f"To convert feet to meters, we multiply by {feet_to_meters_conversion}.")
print(f"Calculation: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")