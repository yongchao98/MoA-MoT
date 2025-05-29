import math

# Volume of the cylindrical barrel
radius_barrel = 10
height_barrel = 15
volume_barrel = math.pi * radius_barrel**2 * height_barrel

# Volume of the hemispherical bucket
radius_bucket = 10
volume_bucket = (2/3) * math.pi * radius_bucket**3

# Number of trips needed
trips = volume_barrel / volume_bucket
trips_rounded = math.ceil(trips)

print(trips_rounded)