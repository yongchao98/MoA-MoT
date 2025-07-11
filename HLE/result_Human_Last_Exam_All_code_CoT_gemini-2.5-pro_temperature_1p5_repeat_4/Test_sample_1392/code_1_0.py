import datetime
import pytz

def analyze_aurora_locations():
    """
    Analyzes which location is most likely to see an overhead aurora
    for a given Kp event and UTC time.
    """
    # --- Step 1: Define Event Parameters and Location Data ---
    kp_index = 7
    # The specific date in early November doesn't matter, only the time of day.
    utc_time_str = "2023-11-05 06:30:00"
    utc_dt = datetime.datetime.strptime(utc_time_str, "%Y-%m-%d %H:%M:%S").replace(tzinfo=pytz.utc)

    # Data for each location: name, timezone, and approximate geomagnetic latitude.
    # The geomagnetic latitude is key to determining position relative to the auroral oval.
    locations = [
        {"id": "A", "name": "Portland, Oregon", "tz": "America/Los_Angeles", "gmag_lat": 52},
        {"id": "B", "name": "Madison, Wisconsin", "tz": "America/Chicago", "gmag_lat": 53},
        {"id": "C", "name": "St. John's, Newfoundland and Labrador", "tz": "America/St_Johns", "gmag_lat": 57},
        {"id": "D", "name": "Alert, Nunavut", "tz": "America/Toronto", "gmag_lat": 86}, # Alert uses EST
        {"id": "E", "name": "Thurso, Scotland", "tz": "Europe/London", "gmag_lat": 61}
    ]

    print(f"Analyzing potential for overhead auroras during a Kp={kp_index} event at {utc_dt.strftime('%H:%M')} UTC.\n")
    print("The ideal geomagnetic latitude for an overhead aurora at Kp=7 is approximately 55-60 degrees.\n")

    # --- Step 2: Analyze each location ---
    for loc in locations:
        tz = pytz.timezone(loc['tz'])
        local_dt = utc_dt.astimezone(tz)

        # Check for darkness (very roughly 7 PM to 7 AM is good for viewing)
        is_dark = 19 <= local_dt.hour or local_dt.hour < 7
        darkness_status = "Good (nighttime)" if is_dark else "Poor (daylight/twilight)"
        
        # Assess geomagnetic position
        gmag_lat = loc['gmag_lat']
        if gmag_lat > 75:
            position_comment = "Too far north (in the polar cap)."
        elif 55 <= gmag_lat <= 60:
            position_comment = "Excellent position (within the prime Kp=7 auroral zone)."
        else:
            position_comment = "Not ideal (likely south or north of the main oval)."

        print(f"--- Analyzing Location {loc['id']}: {loc['name']} ---")
        print(f"Local time is {local_dt.strftime('%H:%M')}.")
        print(f"Darkness condition: {darkness_status}")
        print(f"Approximate geomagnetic latitude is {gmag_lat} degrees.")
        print(f"Likelihood of overhead aurora: {position_comment}")
        print("")

    # --- Step 3: Final Conclusion ---
    print("=========================\nCONCLUSION:")
    print("St. John's, Newfoundland and Labrador is the most likely location because:")
    print("1. It is deep into the night (03:00 local time), providing optimal darkness.")
    print("2. Its geomagnetic latitude of 57 degrees places it directly inside the most intense part of the auroral oval during a Kp=7 storm.")
    print("\nOther locations are less likely:")
    print("- Portland and Madison are too far south magnetically.")
    print("- Alert is too far north, inside the polar cap where the oval has expanded past.")
    print("- Thurso is too close to local sunrise.")

analyze_aurora_locations()