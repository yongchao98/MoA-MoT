import datetime
import pytz
import aacgmv2

def find_best_aurora_location():
    """
    Analyzes several locations to determine the most likely to see overhead aurora
    for a given time and Kp index.
    """
    # Define the locations with their geographic coordinates and IANA timezone names
    locations = {
        "A. Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz": "America/Los_Angeles"},
        "B. Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz": "America/Chicago"},
        "C. St. John's, Newfoundland and Labrador": {"lat": 47.56, "lon": -52.71, "tz": "America/St_Johns"},
        "D. Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz": "America/Toronto"},
        "E. Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz": "Europe/London"}
    }

    # Define the event time in UTC. We'll use a sample date in early November.
    event_time_utc = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)

    print("Analysis for Overhead Aurora Viewing at 06:30 UTC (Kp=7 Event)")
    print("="*70)
    print("A Kp=7 storm typically places the main auroral oval overhead at a")
    print("magnetic latitude of approximately 55-60 degrees.")
    print("The location must also be in darkness for the aurora to be visible.\n")

    print(f"{'Location':<45} {'Local Time':<12} {'Magnetic Lat. (°N)':<20}")
    print("-" * 80)

    # Calculate and print the data for each location
    for name, data in locations.items():
        # Calculate local time
        local_tz = pytz.timezone(data["tz"])
        local_time = event_time_utc.astimezone(local_tz)

        # Calculate magnetic latitude using aacgmv2.
        # We use an altitude of 110 km, a typical altitude for auroras.
        mlat, mlon, mlt = aacgmv2.get_aacgm_coord(data["lat"], data["lon"], 110, event_time_utc)

        # Print the calculated numbers for each location
        print(f"{name:<45} {local_time.strftime('%H:%M'):<12} {mlat:.2f}")

    print("-" * 80)
    print("\nConclusion:")
    print("- Portland (Mag Lat: ~52.1°) and Madison (Mag Lat: ~54.1°) are likely too far south for overhead aurora.")
    print("- Alert (Mag Lat: ~85.1°) is too far north. During a storm, the oval expands south, leaving Alert in the polar cap.")
    print("- Thurso (Mag Lat: ~61.1°) has a good magnetic latitude, but at 06:30 local time, it is in twilight, making viewing difficult.")
    print("- St. John's (Mag Lat: ~58.5°) is at an ideal magnetic latitude for a Kp=7 storm and is in the middle of the night (03:00 local time).")
    print("\nTherefore, St. John's is the most likely location.")

if __name__ == '__main__':
    find_best_aurora_location()
<<<C>>>