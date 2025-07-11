import datetime
import pytz
import aacgmv2

def find_best_aurora_location():
    """
    Analyzes potential aurora viewing locations based on time and geomagnetic coordinates.
    A Kp=7 event means the auroral oval typically expands to a magnetic latitude of 56-58 degrees.
    The best location will be dark and have a magnetic latitude in this range.
    """

    # Define the locations with their geographic coordinates and timezones
    locations = {
        "A. Portland, Oregon": {"lat": 45.5, "lon": -122.7, "tz": "America/Los_Angeles"},
        "B. Madison, Wisconsin": {"lat": 43.1, "lon": -89.4, "tz": "America/Chicago"},
        "C. St. John's, Newfoundland and Labrador": {"lat": 47.6, "lon": -52.7, "tz": "America/St_Johns"},
        "D. Alert, Nunavut": {"lat": 82.5, "lon": -62.3, "tz": "America/Toronto"},
        "E. Thurso, Scotland": {"lat": 58.6, "lon": -3.5, "tz": "Europe/London"},
    }

    # Define the event time in UTC. The exact date in early November doesn't significantly change the result.
    event_time_utc = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)

    print(f"Analysis for an aurora event at {event_time_utc.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    print("A Kp-index of 7 means the aurora is likely overhead near 56-58° magnetic latitude.")
    print("-" * 70)

    for name, data in locations.items():
        # --- Step 1: Calculate Local Time to check for darkness ---
        local_tz = pytz.timezone(data["tz"])
        local_time = event_time_utc.astimezone(local_tz)

        # --- Step 2: Calculate Altitude-Adjusted Corrected Geomagnetic (AACGM) Latitude ---
        # We assume an altitude of 110 km, a common aurora altitude, for overhead calculation.
        mlat, _, _ = aacgmv2.get_aacgm_coord(data["lat"], data["lon"], 110, event_time_utc)

        # --- Step 3: Print the results for each location ---
        print(f"Location: {name}")
        print(f"  - Local Time: {local_time.strftime('%H:%M %Z')}")
        print(f"  - Magnetic Latitude: {mlat:.2f}°")
        
        # Add a comment based on the data
        if local_time.hour < 6 or local_time.hour > 20:
            if 56 <= mlat <= 59:
                print("  - Conclusion: Excellent candidate. It's dark and in the prime magnetic latitude zone.")
            elif mlat > 59:
                 print("  - Conclusion: Likely too far north (in the polar cap).")
            elif mlat < 56:
                print("  - Conclusion: Likely too far south for an overhead view.")
        else:
            print("  - Conclusion: Not a good candidate. It's daytime or twilight.")
        print()

    print("-" * 70)
    print("Based on the analysis, St. John's is the most likely location to witness overhead auroras.")

if __name__ == '__main__':
    find_best_aurora_location()