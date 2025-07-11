# The elevation of the unincorporated community, Descanso, CA, according to the GNIS.
elevation_in_feet = 3465

# The conversion factor from feet to meters.
feet_to_meters_conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor

# Print the full equation for the user.
print(f"The elevation of Descanso, California is {elevation_in_feet} feet.")
print(f"To find the elevation in meters, we multiply by the conversion factor {feet_to_meters_conversion_factor}.")
print(f"Final equation: {elevation_in_feet} * {feet_to_meters_conversion_factor} = {elevation_in_meters}")