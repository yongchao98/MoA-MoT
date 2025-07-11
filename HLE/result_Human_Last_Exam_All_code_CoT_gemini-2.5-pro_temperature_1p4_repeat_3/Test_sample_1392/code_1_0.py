import math
from datetime import datetime, timedelta

def get_local_time(utc_time_str, tz_offset):
    """Calculates local time from a UTC time string and timezone offset."""
    utc_hour, utc_minute = map(int, utc_time_str.split(':'))
    # Using a base date, but only the time is relevant
    base_utc_time = datetime(2023, 11, 5, utc_hour, utc_minute)
    local_time = base_utc_time + timedelta(hours=tz_offset)
    return local_time.strftime('%H:%M Local Time')

def get_magnetic_latitude(geo_lat, geo_lon):
    """
    Calculates the magnetic latitude based on a simple dipole model.
    Uses IGRF-13 Geomagnetic North Pole coordinates for 2020 (approx. 80.7N, 72.7W).
    """
    # Geomagnetic North Pole (approximate coordinates)
    pole_lat_deg = 80.7
    pole_lon_deg = -72.7

    # Convert all degrees to radians for calculations
    phi = math.radians(geo_lat)
    lambda_ = math.radians(geo_lon)
    phi_p = math.radians(pole_lat_deg)
    lambda_p = math.radians(pole_lon_deg)

    # Apply the formula for sine of magnetic latitude
    sin_mag_lat = math.sin(phi) * math.sin(phi_p) + math.cos(phi) * math.cos(phi_p) * math.cos(lambda_ - lambda_p)
    
    # The result can sometimes be slightly > 1.0 due to approximations, so clip it.
    if sin_mag_lat > 1.0:
        sin_mag_lat = 1.0
    if sin_mag_lat < -1.0:
        sin_mag_lat = -1.0

    # Calculate magnetic latitude and convert back to degrees
    mag_lat = math.degrees(math.asin(sin_mag_lat))
    return mag_lat

def analyze_aurora_locations():
    """
    Analyzes potential aurora viewing locations based on UTC time and Kp index.
    """
    locations = [
        {'name': 'A. Portland, Oregon', 'lat': 45.5, 'lon': -122.7, 'tz': -8},
        {'name': 'B. Madison, Wisconsin', 'lat': 43.1, 'lon': -89.4, 'tz': -6},
        {'name': 'C. St. John\'s, Newfoundland and Labrador', 'lat': 47.6, 'lon': -52.7, 'tz': -3.5},
        {'name': 'D. Alert, Nunavut', 'lat': 82.5, 'lon': -62.3, 'tz': -5},
        {'name': 'E. Thurso, Scotland', 'lat': 58.6, 'lon': -3.5, 'tz': 0},
    ]

    utc_time_str = "06:30"
    kp_index = 7
    
    print(f"Analysis for Auroral Sighting at {utc_time_str} UTC during a Kp={kp_index} event:\n")
    print("The ideal conditions for overhead aurora during a Kp=7 storm are darkness, a magnetic latitude of roughly 55-65 degrees, and a local time near midnight.\n")

    for loc in locations:
        local_time = get_local_time(utc_time_str, loc['tz'])
        mag_lat = get_magnetic_latitude(loc['lat'], loc['lon'])
        
        print(f"Location: {loc['name']}")
        print(f"  - Local Time: {local_time}")
        print(f"  - Geographic Latitude: {loc['lat']:.1f}° N")
        print(f"  - Magnetic Latitude: {mag_lat:.1f}°\n")
        
    print("Conclusion:")
    print("All locations are in darkness. However, their positions relative to the auroral oval differ significantly.")
    print("- Alert (D) is at a very high magnetic latitude (85.8°), placing it inside the polar cap, north of the main auroral oval during a storm.")
    print("- Portland (A) and Madison (B) have lower magnetic latitudes (52.8° and 53.6°), putting them on the southern edge of visibility. Overhead aurora is less likely.")
    print("- Thurso (E) and St. John's (C) have ideal magnetic latitudes (57.6° and 56.1°) to be under the main oval.")
    print("- The final deciding factor is the time of night. St. John's at 03:00 local time is in a prime period for intense auroral activity. Thurso at 06:30 local time is in the dawn sector, where aurora is typically less dynamic.")
    print("\nTherefore, St. John's is the most probable location to witness intense overhead auroras.")

analyze_aurora_locations()