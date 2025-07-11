# Step 1: Define the elevation in feet from the Geographic Names Information System (GNIS).
# The elevation for Jacumba, California is 2805 feet.
elevation_feet = 2805

# Step 2: Define the conversion factor from feet to meters.
feet_to_meters_conversion = 0.3048

# Step 3: Calculate the elevation in meters.
elevation_meters = elevation_feet * feet_to_meters_conversion

# Step 4: Print the full equation and the final result.
print(f"The elevation of the community is {elevation_feet} feet.")
print(f"The conversion factor is {feet_to_meters_conversion} meters per foot.")
print(f"The elevation in meters is calculated as: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")