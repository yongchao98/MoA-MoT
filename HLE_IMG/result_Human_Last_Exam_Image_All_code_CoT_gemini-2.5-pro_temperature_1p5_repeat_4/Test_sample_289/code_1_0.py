# The elevation of Boulevard, California, according to the GNIS is 3,668 feet.
# This script converts the elevation from feet to meters.

# Elevation in feet
elevation_feet = 3668

# Conversion factor from feet to meters
feet_to_meters_conversion = 0.3048

# Calculate the elevation in meters
elevation_meters = elevation_feet * feet_to_meters_conversion

# Print the full equation and the result
print(f"The official elevation of the unincorporated community is {elevation_feet} feet.")
print("To convert feet to meters, we multiply by the conversion factor.")
print(f"Calculation: {elevation_feet} feet * {feet_to_meters_conversion} meters/foot = {elevation_meters:.2f} meters")

# Round the final answer to the nearest whole number for the final output
rounded_elevation_meters = round(elevation_meters)
print(f"The elevation is approximately {rounded_elevation_meters} meters.")