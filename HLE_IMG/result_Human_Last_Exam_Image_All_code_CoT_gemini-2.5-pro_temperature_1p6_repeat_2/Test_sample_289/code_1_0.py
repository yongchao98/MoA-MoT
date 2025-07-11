# The elevation of Descanso, CA according to GNIS data is 3465 feet.
elevation_in_feet = 3465

# The conversion factor from feet to meters.
feet_to_meters_conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor

print(f"The elevation of the unincorporated community is {elevation_in_feet} feet.")
print(f"To convert feet to meters, we multiply by {feet_to_meters_conversion_factor}.")
print(f"The final calculation is: {elevation_in_feet} * {feet_to_meters_conversion_factor} = {elevation_in_meters}")
print(f"The elevation in meters is: {elevation_in_meters}")