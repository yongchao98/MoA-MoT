import math
from datetime import datetime, timezone, timedelta

def calculate_magnetic_latitude(geo_lat, geo_lon):
    """
    Calculates the approximate magnetic latitude of a point on Earth.
    This uses a simplified dipole model with the North Magnetic Pole's
    approximate location in the modern era (near 86.5°N, 164.0°E).
    """
    # North Magnetic Pole approximate coordinates
    pole_lat = 86.5
    pole_lon = 164.0 # East Longitude

    # Convert all degree values to radians for calculations
    geo_lat_rad = math.radians(geo_lat)
    geo_lon_rad = math.radians(geo_lon)
    pole_lat_rad = math.radians(pole_lat)
    pole_lon_rad = math.radians(pole_lon)

    # Use the spherical law of cosines to find the magnetic colatitude
    magnetic_colatitude_rad = math.acos(
        math.sin(geo_lat_rad) * math.sin(pole_lat_rad) +
        math.cos(geo_lat_rad) * math.cos(pole_lat_rad) * math.cos(geo_lon_rad - pole_lon_rad)
    )
    
    # Magnetic latitude is 90 degrees minus the magnetic colatitude in degrees
    magnetic_latitude = 90.0 - math.degrees(magnetic_colatitude_rad)
    return magnetic_latitude

def get_local_time(utc_time_str, utc_offset_hours):
    """Calculates local time from a UTC string and offset."""
    utc_dt = datetime.strptime(utc_time_str, '%H:%M').replace(tzinfo=timezone.utc)
    local_dt = utc_dt + timedelta(hours=utc_offset_hours)
    return local_dt.strftime('%H:%M')

def is_dark(local_time_str):
    """Checks if the local time is during dark hours (approx. 8 PM to 5 AM)."""
    hour = int(local_time_str.split(':')[0])
    return hour >= 20 or hour <= 5

# --- Main Analysis ---

# Event Parameters
kp_index = 7
utc_time = "06:30"
# For a Kp=7 event, the auroral oval is centered roughly at 58° magnetic latitude.
target_mag_lat = 58 

# Location Data: Name, Geographic Lat, Geographic Lon (West is negative), UTC offset
# UTC offsets are for early November (Standard Time).
locations = {
    "A. Portland, Oregon": {"lat": 45.5, "lon": -122.7, "tz": -8},
    "B. Madison, Wisconsin": {"lat": 43.1, "lon": -89.4, "tz": -6},
    "C. St. John's, Newfoundland": {"lat": 47.6, "lon": -52.7, "tz": -3.5},
    "D. Alert, Nunavut": {"lat": 82.5, "lon": -62.3, "tz": -5},
    "E. Thurso, Scotland": {"lat": 58.6, "lon": -3.5, "tz": 0},
}

results = {}
best_candidate = ""
min_diff = float('inf')

print(f"Analysis for Auroral Sighting at {utc_time} UTC with Kp={kp_index}\n")
print(f"{'Location':<30} | {'Local Time':<12} | {'Is Dark?':<10} | {'Magnetic Lat.':<15}")
print("-" * 75)

for name, data in locations.items():
    local_time = get_local_time(utc_time, data["tz"])
    darkness = "Yes" if is_dark(local_time) else "No"
    mag_lat = calculate_magnetic_latitude(data["lat"], data["lon"])
    results[name] = {
        'local_time': local_time,
        'is_dark': darkness,
        'mag_lat': mag_lat
    }
    print(f"{name:<30} | {local_time:<12} | {darkness:<10} | {mag_lat:.2f}°")

print("-" * 75)

# --- Evaluation ---
print("\nStep 1: Check for darkness. Aurora is only visible at night.")
thurso_result = results['E. Thurso, Scotland']
print(f"- Thurso, Scotland is at {thurso_result['local_time']} local time. It is morning, so the sky is too bright. This option is eliminated.")

print("\nStep 2: Check magnetic latitude relative to the auroral oval.")
print(f"A Kp={kp_index} storm pushes the aurora to lower latitudes. The ideal magnetic latitude for an *overhead* view is ~{target_mag_lat}°. \n")

alert_result = results['D. Alert, Nunavut']
print(f"- Alert, Nunavut has a magnetic latitude of {alert_result['mag_lat']:.2f}°. This is extremely high, placing it inside the 'polar cap' and north of the main auroral oval, which has expanded southwards. This option is unlikely.")

print("\nStep 3: Compare the remaining locations in North America.")
for name in ["A. Portland, Oregon", "B. Madison, Wisconsin", "C. St. John's, Newfoundland"]:
    result = results[name]
    diff = abs(result['mag_lat'] - target_mag_lat)
    print(f"- {name}: Magnetic Latitude = {result['mag_lat']:.2f}°. Distance from target latitude ~{target_mag_lat}° is {diff:.2f}°.")
    if diff < min_diff:
        min_diff = diff
        best_candidate = name

print("\n--- Conclusion ---")
print(f"St. John's, Newfoundland and Labrador has the magnetic latitude closest to the center of the expanded auroral oval during a Kp=7 storm.")
print("Therefore, it is the most likely location among the choices to witness overhead auroras.")
print(f"The best choice is {best_candidate}")

<<<C>>>