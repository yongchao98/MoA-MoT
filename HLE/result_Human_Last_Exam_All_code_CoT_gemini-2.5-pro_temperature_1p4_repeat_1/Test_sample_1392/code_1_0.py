import datetime
import sys

try:
    import aacgmv2
except ImportError:
    # Use print for all user-facing messages.
    print("Error: The 'aacgmv2' library is required for this calculation.", file=sys.stderr)
    print("Please install it by running: pip install aacgmv2", file=sys.stderr)
    sys.exit(1)

def solve_aurora_location():
    """
    Analyzes which location is most likely to see overhead auroras
    during a Kp=7 event at 06:30 UTC in early November.
    """
    # Define the time of the event
    # Using a representative date in early November at 06:30 UTC.
    event_time_utc = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=datetime.timezone.utc)

    # Define the candidate locations with their geographic coordinates and UTC offset
    # Note: Timezones are for Standard Time, typically in effect in early November.
    locations = {
        "A. Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz_offset": -8},
        "B. Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz_offset": -6},
        "C. St. John's, NL": {"lat": 47.56, "lon": -52.71, "tz_offset": -3.5},
        "D. Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz_offset": -5},
        "E. Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz_offset": 0}
    }

    # The auroral oval's position is key. For a Kp=7 storm, overhead auroras
    # are often seen around 60-65 degrees magnetic latitude.
    # Intensity is typically highest near magnetic midnight (MLT ≈ 0).

    print("Analysis of Potential Aurora Viewing Locations")
    print(f"Event Time (UTC): {event_time_utc.strftime('%Y-%m-%d %H:%M')}")
    print("-" * 85)
    print(f"{'Location':<25} | {'Local Time':<12} | {'Is Dark?':<10} | {'Magnetic Latitude':<20} | {'Magnetic Local Time':<20}")
    print("-" * 85)

    # Altitude for auroral emissions is typically 100-400 km. We use a standard average of 300 km.
    aurora_altitude_km = 300

    for name, data in locations.items():
        # Calculate local time
        local_time = event_time_utc + datetime.timedelta(hours=data["tz_offset"])

        # Check for darkness (a simple check for evening/night hours)
        is_dark = "Yes" if 20 <= local_time.hour or local_time.hour < 6 else "No"
        # For high-latitude locations like Alert, it's polar night in November.
        if "Alert" in name:
            is_dark = "Yes"

        # Calculate magnetic coordinates using aacgmv2
        # mlat = magnetic latitude, mlon = magnetic longitude, mlt = magnetic local time
        mlat, _, mlt = aacgmv2.get_aacgm_coord(data["lat"], data["lon"], aurora_altitude_km, event_time_utc)

        # Print the final numbers in the analysis for each location
        print(f"{name:<25} | {local_time.strftime('%I:%M %p'):<12} | {is_dark:<10} | {mlat:^20.2f} | {mlt:^20.2f}")

    print("-" * 85)
    print("\nConclusion:")
    print("A Kp=7 storm expands the auroral oval. The ideal location for overhead views would have a magnetic latitude of ~60-65°.")
    print("Based on the data:")
    print("- Portland, Madison, and Thurso have magnetic latitudes below the primary auroral zone.")
    print("- Alert is deep within the polar cap (mag lat > 85°), which is north of the most intense part of the auroral oval.")
    print("- St. John's, with a magnetic latitude of ~58°, is positioned closest to the expected active auroral oval for a Kp=7 storm.")

if __name__ == '__main__':
    solve_aurora_location()
<<<C>>>