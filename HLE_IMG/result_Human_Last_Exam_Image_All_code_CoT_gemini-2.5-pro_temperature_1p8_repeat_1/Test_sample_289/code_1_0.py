# The elevation of the unincorporated community, Warner Springs, CA,
# is 3136 feet according to the Geographic Names Information System (GNIS).
# This script converts that elevation to meters.

# The elevation in feet.
elevation_feet = 3136

# The conversion factor from feet to meters.
conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_meters = elevation_feet * conversion_factor

# Print the context and the full calculation.
print(f"The elevation of Warner Springs, CA is {elevation_feet} feet.")
print(f"The conversion calculation to meters is as follows:")
print(f"{elevation_feet} * {conversion_factor} = {elevation_meters}")
print(f"The elevation is {elevation_meters} meters.")