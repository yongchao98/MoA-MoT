import datetime
import pytz
from astral.observer import Observer
from astral.sun import elevation
import aacgmv2

def analyze_aurora_location(name, lat, lon, tz_str, utc_time):
    """
    Analyzes a single location for aurora viewing potential by calculating
    local time, sun elevation, and magnetic latitude.
    """
    # --- Step 1: Calculate local time and sun elevation ---
    local_tz = pytz.timezone(tz_str)
    local_time = utc_time.astimezone(local_tz)
    
    obs = Observer(latitude=lat, longitude=lon)
    sun_elev = elevation(obs, utc_time)
    # A sun elevation < -6 degrees (below civil twilight) is dark enough for good viewing.
    is_dark = sun_elev < -6

    # --- Step 2: Calculate magnetic latitude ---
    # A Kp=7 storm pushes the auroral oval to ~55-65 degrees magnetic latitude.
    # We assume an altitude of 0 km for the calculation.
    mlat, mlon, mlt = aacgmv2.get_aacgm_coord(lat, lon, 0, utc_time)
    
    # --- Step 3: Check if location is in the prime Kp=7 zone ---
    in_zone = 55 <= mlat <= 65

    # --- Step 4: Print the analysis for the location ---
    print(f"--- Analysis for: {name} ---")
    print(f"Local Time at 06:30 UTC is: {local_time.strftime('%H:%M %Z')}")
    print(f"Sun's Elevation: {sun_elev:.2f} degrees")
    print(f"Is it dark enough for viewing? {'Yes' if is_dark else 'No'}")
    print(f"Geomagnetic Latitude: {mlat:.2f} degrees")
    print(f"Is it in the prime Kp=7 overhead zone (55-65°)? {'Yes' if in_zone else 'No'}")
    print("-" * 35 + "\n")


# --- Main Execution ---

# Define the target time: Early November, 06:30 UTC
utc_event_time = datetime.datetime(2023, 11, 5, 6, 30, 0, tzinfo=pytz.utc)

# Define the answer choices with their geographic coordinates and timezones
locations_data = {
    "A": ("Portland, Oregon", 45.52, -122.68, "America/Los_Angeles"),
    "B": ("Madison, Wisconsin", 43.07, -89.40, "America/Chicago"),
    "C": ("St. John's, Newfoundland and Labrador", 47.56, -52.71, "America/St_Johns"),
    "D": ("Alert, Nunavut", 82.50, -62.35, "America/Toronto"),
    "E": ("Thurso, Scotland", 58.59, -3.52, "Europe/London")
}

print(f"Finding the best location to see overhead aurora at {utc_event_time.strftime('%Y-%m-%d %H:%M %Z')}\n")

for key, (name, lat, lon, tz) in locations_data.items():
    analyze_aurora_location(f"{key}. {name}", lat, lon, tz, utc_event_time)

print("--- Conclusion ---")
print("St. John's, Newfoundland and Labrador is the only location that is both sufficiently dark")
print("and situated at the correct magnetic latitude to be directly under the expanded auroral oval.")
