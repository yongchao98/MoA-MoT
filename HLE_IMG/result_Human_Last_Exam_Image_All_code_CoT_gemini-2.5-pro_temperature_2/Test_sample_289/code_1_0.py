# The elevation of Descanso, CA is 3514 feet according to the GNIS.
elevation_feet = 3514

# The conversion factor from feet to meters.
feet_to_meters_conversion = 0.3048

# Calculate the elevation in meters.
elevation_meters = elevation_feet * feet_to_meters_conversion

# Print the equation with all the numbers.
print(f"The elevation of the unincorporated community is {elevation_feet} feet.")
print(f"To convert to meters, we use the formula: Elevation in feet * {feet_to_meters_conversion}.")
print(f"Calculation: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")
print(f"The elevation is approximately {round(elevation_meters, 2)} meters.")