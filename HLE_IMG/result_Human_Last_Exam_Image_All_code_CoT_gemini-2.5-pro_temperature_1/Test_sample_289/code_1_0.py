# Step 1: Define the elevation of Descanso, CA in feet, as found in the GNIS database.
elevation_feet = 3527

# Step 2: Define the conversion factor from feet to meters.
feet_to_meters_conversion = 0.3048

# Step 3: Calculate the elevation in meters.
elevation_meters = elevation_feet * feet_to_meters_conversion

# Step 4: Print the equation and the final result.
print(f"The elevation of Descanso, CA is {elevation_feet} feet.")
print("To convert feet to meters, we use the conversion factor: 1 foot = 0.3048 meters.")
print("Calculation: {} feet * {} meters/foot = {} meters".format(elevation_feet, feet_to_meters_conversion, elevation_meters))
print(f"The elevation of Descanso, CA is approximately {round(elevation_meters)} meters.")