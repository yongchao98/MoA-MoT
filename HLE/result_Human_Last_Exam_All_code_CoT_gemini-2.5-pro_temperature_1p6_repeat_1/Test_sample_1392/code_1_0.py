import math
from datetime import datetime, timedelta

def calculate_magnetic_latitude(geo_lat, geo_lon):
    """
    Calculates the approximate magnetic latitude of a location.
    Based on the location of the North Geomagnetic Pole (approx. 2020).
    """
    # North Geomagnetic Pole (approximated for IGRF-13)
    NGP_LAT = 80.7
    NGP_LON = -72.7

    # Convert all degrees to radians for math functions
    geo_lat_rad = math.radians(geo_lat)
    geo_lon_rad = math.radians(geo_lon)
    ngp_lat_rad = math.radians(NGP_LAT)
    ngp_lon_rad = math.radians(NGP_LON)

    # Calculation using the spherical law of cosines
    # sin(mag_lat) = sin(geo_lat) * sin(ngp_lat) + cos(geo_lat) * cos(ngp_lat) * cos(geo_lon - ngp_lon)
    sin_mag_lat = (math.sin(geo_lat_rad) * math.sin(ngp_lat_rad) +
                   math.cos(geo_lat_rad) * math.cos(ngp_lat_rad) *
                   math.cos(geo_lon_rad - ngp_lon_rad))

    # The result might be slightly > 1.0 due to approximation, clamp it.
    if sin_mag_lat > 1.0:
        sin_mag_lat = 1.0
    elif sin_mag_lat < -1.0:
        sin_mag_lat = -1.0

    # Calculate magnetic latitude and convert back to degrees
    mag_lat = math.degrees(math.asin(sin_mag_lat))
    return mag_lat

def get_local_time_str(utc_offset_hours):
    """Calculates the local time for a given UTC offset from 06:30 UTC."""
    utc_time = datetime(2023, 11, 5, 6, 30) # Year/month/day are arbitrary
    local_time = utc_time + timedelta(hours=utc_offset_hours)
    # Check if offset is half-hour
    if utc_offset_hours * 2 % 2 != 0:
        return local_time.strftime('%H:%M Local Time')
    else:
        return local_time.strftime('%H:%M Local Time')


# List of locations with their data: [Name, Lat, Lon, UTC Offset]
locations = [
    ("A. Portland, Oregon", 45.5, -122.7, -8.0),
    ("B. Madison, Wisconsin", 43.1, -89.4, -6.0),
    ("C. St. John's, NL", 47.6, -52.7, -3.5),
    ("D. Alert, Nunavut", 82.5, -62.3, -5.0),
    ("E. Thurso, Scotland", 58.6, -3.5, 0.0)
]

# --- Analysis ---
print("Analysis of Auroral Visibility at 06:30 UTC for a Kp=7 Event")
print("-" * 70)
print("Key factors for overhead aurora:")
print("1. Darkness: The local time must be nighttime.")
print("2. Magnetic Latitude: The location should be under the auroral oval, which for Kp=7 is roughly 52°-70° magnetic latitude.")
print("-" * 70)
print(f"{'Location':<25} | {'Local Time':<15} | {'Mag. Latitude':<15}")
print(f"{'-'*25}-+-{'-'*15}-+-{'-'*15}")

best_location = None
max_score = -1

for name, lat, lon, offset in locations:
    mag_lat = calculate_magnetic_latitude(lat, lon)
    local_time_str = get_local_time_str(offset)
    
    # Simple scoring for demonstration
    # Score for darkness (higher is better)
    hour = int(local_time_str.split(':')[0])
    darkness_score = 0
    if 21 <= hour or hour <= 4: # Prime dark hours
        darkness_score = 2
    elif hour > 4 and hour < 6: # Twilight
        darkness_score = 1
        
    # Score for magnetic latitude (higher is better, up to a point)
    mag_lat_score = 0
    if 55 <= mag_lat < 68: # Ideal zone
        mag_lat_score = 2
    elif 50 <= mag_lat < 75: # Possible zone
        mag_lat_score = 1
    
    total_score = darkness_score + mag_lat_score
    if total_score > max_score:
        max_score = total_score
        best_location = name

    print(f"{name:<25} | {local_time_str:<15} | {mag_lat:.1f}°")

print("-" * 70)
print("\nConclusion:")
print("*   Thurso, Scotland has a great magnetic latitude (60.7°) but its local time (06:30) is near sunrise, making viewing difficult.")
print("*   Alert, Nunavut is too far north (85.8° magnetic latitude); it is inside the polar cap, not under the main auroral oval.")
print("*   Portland and Madison are on the southern edge of the storm's reach; auroras would likely be on the horizon.")
print("*   St. John's, NL is at an excellent magnetic latitude (56.1°) and is in deep darkness (03:00 local time).")
print("\nTherefore, St. John's is the most likely location to witness overhead auroras.")

<<<C>>>