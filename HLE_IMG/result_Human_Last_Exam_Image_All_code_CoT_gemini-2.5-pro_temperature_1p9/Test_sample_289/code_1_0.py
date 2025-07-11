# The elevation of Potrero, CA is 2319 feet according to the GNIS.
elevation_feet = 2319

# The conversion factor from feet to meters is 0.3048.
feet_to_meters_conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_meters = elevation_feet * feet_to_meters_conversion_factor

# Print the equation and the result.
print(f"The elevation of the unincorporated community in feet is: {elevation_feet}")
print(f"The conversion factor from feet to meters is: {feet_to_meters_conversion_factor}")
print(f"The calculation is: {elevation_feet} feet * {feet_to_meters_conversion_factor} meters/foot")
print(f"The elevation in meters is: {elevation_meters}")