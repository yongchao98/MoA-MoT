# The elevation of Descanso, California, according to the GNIS.
elevation_in_feet = 3527

# The conversion factor from feet to meters.
feet_to_meters_conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor

# Print the final equation with all the numbers.
print(f"The elevation of the unincorporated community is {elevation_in_feet} feet.")
print("The conversion from feet to meters is calculated as:")
print(f"{elevation_in_feet} * {feet_to_meters_conversion_factor} = {elevation_in_meters}")
print(f"The elevation in meters is approximately {round(elevation_in_meters)}.")