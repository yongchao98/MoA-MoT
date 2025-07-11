# First, ensure required libraries are installed:
# pip install aacgmv2==2.6.2 pytz==2023.3

import datetime
import pytz
import aacgmv2

def find_best_aurora_location():
    """
    Analyzes which location is most likely to see overhead auroras
    for a given Kp event and time.
    """
    # Define event parameters
    kp_index = 7
    # Using a sample date in early November. The result is not sensitive to the exact day.
    event_utc_time = datetime.datetime(2023, 11, 5, 6, 30, 0, tzinfo=pytz.utc)

    # Define the answer choices with their geographic coordinates and timezones
    locations = {
        "A. Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz": "America/Los_Angeles"},
        "B. Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz": "America/Chicago"},
        "C. St. John's, Newfoundland and Labrador": {"lat": 47.56, "lon": -52.71, "tz": "America/St_Johns"},
        "D. Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz": "America/Toronto"},
        "E. Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz": "Europe/London"},
    }

    # For a given Kp index, the center of the auroral oval is roughly at this magnetic latitude.
    # Formula: Magnetic Latitude ≈ 66.5 - Kp
    target_mlat = 66.5 - kp_index

    print(f"Analysis for Kp={kp_index} event at {event_utc_time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    print(f"The ideal magnetic latitude for overhead aurora is approximately {target_mlat:.1f}°.\n")

    results = []
    for name, data in locations.items():
        # Calculate magnetic latitude (mlat) and magnetic local time (mlt) for the event
        # 300km is a typical altitude for auroral displays.
        mlat, mlon, mlt = aacgmv2.get_aacgm_coord(data["lat"], data["lon"], 300, event_utc_time)

        # Calculate local time for context
        local_tz = pytz.timezone(data["tz"])
        local_time = event_utc_time.astimezone(local_tz)

        # Store results for comparison
        results.append({
            "name": name,
            "mlat": mlat,
            "mlt": mlt,
            "local_time_str": local_time.strftime('%H:%M %Z'),
            "diff_from_target": abs(mlat - target_mlat)
        })

    # Find the location with the minimum difference from the target latitude
    best_location = min(results, key=lambda x: x['diff_from_target'])

    # Print the analysis for each location
    print("--- Location Analysis ---")
    for res in sorted(results, key=lambda x: x['diff_from_target']):
        print(f"Location: {res['name']}")
        print(f"  - Local Time: {res['local_time_str']}")
        print(f"  - Magnetic Latitude: {res['mlat']:.1f}° (Difference from target: {res['diff_from_target']:.1f}°)")
        print(f"  - Magnetic Local Time (MLT): {res['mlt']:.1f} (Peak activity is near 00:00)\n")

    print("--- Conclusion ---")
    print(f"The location with the magnetic latitude closest to the ideal overhead viewing zone ({target_mlat:.1f}°) is {best_location['name']}.")
    print(f"Its magnetic local time ({best_location['mlt']:.1f}) is also in the prime post-midnight sector for strong auroral displays.")

if __name__ == '__main__':
    find_best_aurora_location()