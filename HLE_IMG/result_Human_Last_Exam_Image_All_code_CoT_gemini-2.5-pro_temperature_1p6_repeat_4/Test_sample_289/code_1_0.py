# The elevation of Santa Ysabel, CA, according to the GNIS is 2986 feet.
elevation_feet = 2986

# The conversion factor from feet to meters is 0.3048.
conversion_factor = 0.3048

# Calculate the elevation in meters.
elevation_meters = elevation_feet * conversion_factor

print(f"The unincorporated community is Santa Ysabel, CA.")
print(f"Its elevation in feet is: {elevation_feet}")
print(f"The conversion factor from feet to meters is: {conversion_factor}")
print(f"The elevation in meters is {elevation_feet} * {conversion_factor} = {elevation_meters:.2f}")
print(f"Rounding to the nearest whole number, the elevation is {round(elevation_meters)} meters.")