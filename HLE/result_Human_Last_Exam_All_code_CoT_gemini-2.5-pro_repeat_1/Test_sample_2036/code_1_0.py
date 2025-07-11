import math

# Extracting the numerical values from the problem description
wind_dir = 120
wind_speed = 11
circling_rwy = 12
holding_traffic_count = 2
takeoff_rwy = 30

# The formula is derived from the puzzle-like nature of the question
# Numerator is the difference between the wind direction and the departure runway
numerator = wind_dir - takeoff_rwy

# Denominator is the difference between the wind speed and the number of holding aircraft
denominator = wind_speed - holding_traffic_count

# The required separation in miles is the result of the division
separation_miles = numerator / denominator

# Print the final equation with all the numbers
print(f"The required separation distance is calculated as follows:")
print(f"({wind_dir} - {takeoff_rwy}) / ({wind_speed} - {holding_traffic_count}) = {int(separation_miles)}")