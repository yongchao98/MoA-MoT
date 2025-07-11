import math

# It is recommended to install the geopy library first
# You can install it by running: pip install geopy
from geopy.distance import geodesic

def dms_to_dd(is_positive, degrees, minutes, seconds):
    """Converts Degrees, Minutes, Seconds to signed Decimal Degrees."""
    dd = degrees + minutes / 60 + seconds / 3600
    if not is_positive:
        dd = -dd
    return dd

# --- Coordinates of the Rock Carving ---
# Provided coordinates: 29° 3' 28.15''N and 103° 48' 11.84''W
rock_lat_d, rock_lat_m, rock_lat_s = 29, 3, 28.15
rock_lon_d, rock_lon_m, rock_lon_s = 103, 48, 11.84

# Convert to Decimal Degrees. North is positive, West is negative.
rock_lat_dd = dms_to_dd(True, rock_lat_d, rock_lat_m, rock_lat_s)
rock_lon_dd = dms_to_dd(False, rock_lon_d, rock_lon_m, rock_lon_s)
rock_coords = (rock_lat_dd, rock_lon_dd)

print("--- Analysis of Geographic Claims ---")
print(f"Rock Carving Latitude ({rock_lat_d}°{rock_lat_m}'{rock_lat_s}''N) = {rock_lat_dd:.4f}°")
print(f"Rock Carving Longitude ({rock_lon_d}°{rock_lon_m}'{rock_lon_s}''W) = {rock_lon_dd:.4f}°\n")


# --- Define locations to check against the claims ---
locations_to_check = {
    "B": {
        "name": "Chisos Mountains",
        "coords": (29.2694, -103.3005),  # Coordinates for a central point in the range
        "claim_dist": 20,
        "claim_dir": "north"
    },
    "D": {
        "name": "Bravo River near Lajitas, TX",
        "coords": (29.2638, -103.7663),  # Coordinates for Lajitas, TX
        "claim_dist": 10,
        "claim_dir": "northwest"
    }
}

# --- Perform and Print Analysis ---
for option, data in locations_to_check.items():
    print(f"--- Verifying Claim in Option {option} ---")
    print(f"Feature: {data['name']}")
    print(f"Claim: Located {data['claim_dist']} miles {data['claim_dir']} from the carving.")

    # Calculate distance
    calculated_dist = geodesic(rock_coords, data['coords']).miles
    
    # Determine direction
    lat_diff = data['coords'][0] - rock_coords[0]
    lon_diff = data['coords'][1] - rock_coords[1]
    
    direction = ""
    if lat_diff > 0:
        direction += "North"
    elif lat_diff < 0:
        direction += "South"
    
    if lon_diff > 0:
        direction += "East"
    elif lon_diff < 0:
        direction += "West"

    # Print results
    print(f"Calculation Result:")
    # Printing each number for the comparison as requested
    print(f"Calculated distance = {calculated_dist:.1f} miles (Claim was {data['claim_dist']} miles)")
    print(f"Calculated direction = {direction} (Claim was {data['claim_dir']})")
    
    is_correct = (abs(calculated_dist - data['claim_dist']) < 5) and (data['claim_dir'].upper() in direction.upper())
    print(f"Conclusion: Claim in Option {option} is {'LIKELY TRUE' if is_correct else 'FALSE'}.\n")
