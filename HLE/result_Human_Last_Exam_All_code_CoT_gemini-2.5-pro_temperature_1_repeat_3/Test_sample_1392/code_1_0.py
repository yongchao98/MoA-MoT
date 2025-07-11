import datetime
import pytz

def solve_aurora_location():
    """
    Analyzes potential aurora viewing locations based on Kp index and UTC time.
    """
    # --- Input Data ---
    # Kp index of the geomagnetic storm
    kp_index = 7
    # UTC time of the event
    utc_time_str = "06:30"

    # Location data: [City, Geographic Latitude, Timezone String, Approx. Magnetic Latitude]
    locations = [
        ["Portland, Oregon", 45.5, "America/Los_Angeles", 51],
        ["Madison, Wisconsin", 43.1, "America/Chicago", 52],
        ["St. John's, Newfoundland and Labrador", 47.6, "America/St_Johns", 56],
        ["Alert, Nunavut", 82.5, "America/Toronto", 86], # Note: Alert uses EST
        ["Thurso, Scotland", 58.6, "Europe/London", 61]
    ]

    # --- Analysis ---
    print(f"Analyzing best location for overhead aurora during a Kp={kp_index} event at {utc_time_str} UTC.\n")

    # The auroral oval's equatorward boundary expands south during a storm.
    # A common estimation for the boundary is: 66.5 - 2 * Kp
    equatorward_boundary = 66.5 - 2 * kp_index
    print(f"For Kp={kp_index}, the estimated equatorward boundary of the auroral oval is ~{equatorward_boundary:.1f}° magnetic latitude.")
    print("Locations need to be dark and have a magnetic latitude at or above this boundary for potential overhead aurora.\n")

    # Parse the UTC time
    utc_hour, utc_minute = map(int, utc_time_str.split(':'))
    # Use a fixed date; the specific day in early November doesn't change the time zone calculations
    event_utc_time = datetime.datetime(2023, 11, 5, utc_hour, utc_minute, tzinfo=pytz.utc)

    best_candidate = None
    best_score = -1

    print("--- Location Analysis ---")
    for name, geo_lat, tz_str, mag_lat in locations:
        # 1. Calculate Local Time
        local_tz = pytz.timezone(tz_str)
        local_time = event_utc_time.astimezone(local_tz)
        
        # 2. Assess Darkness
        # Ideal viewing is between 8 PM (20:00) and 5 AM (05:00)
        is_dark = 20 <= local_time.hour or local_time.hour < 5
        # Special case for Alert in November (polar night)
        if name == "Alert, Nunavut":
            is_dark = True
        
        # 3. Assess Magnetic Position
        # Score based on how well it fits the criteria
        score = 0
        notes = []
        if is_dark:
            score += 1
            notes.append("Dark sky.")
        else:
            notes.append("Daylight/twilight.")

        # Check magnetic latitude
        if mag_lat < equatorward_boundary:
            notes.append(f"South of the auroral oval (Mag Lat: {mag_lat}°).")
        elif mag_lat > 75: # Deep in polar cap
            score += 0.5 # Can have polar cap aurora, but less likely for main oval
            notes.append(f"In the polar cap, likely north of the main oval (Mag Lat: {mag_lat}°).")
        else: # In the prime zone
            score += 1
            notes.append(f"Excellent magnetic position (Mag Lat: {mag_lat}°).")

        if score > best_score:
            best_score = score
            best_candidate = name

        print(f"Location: {name}")
        print(f"  - Local Time: {local_time.strftime('%H:%M %Z')}")
        print(f"  - Analysis: {' '.join(notes)}\n")

    print("--- Conclusion ---")
    print(f"The location most likely to witness *overhead* auroras is St. John's, Newfoundland and Labrador.")
    print("It is in the middle of the night (optimal darkness) and its magnetic latitude of 56° places it squarely under the main part of the auroral oval during a Kp=7 storm.")

solve_aurora_location()
<<<C>>>