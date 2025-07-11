import datetime
import pytz
import aacgmv2
from astral.geocoder import geocoder
from astral.sun import sun

def solve_aurora_location():
    """
    Analyzes which location is most likely to see overhead auroras for a given
    geomagnetic storm level (Kp-index) and time.
    """
    # --- Input Data ---
    # The event details from the user's question
    kp_index = 7
    # We use a representative date in early November. The exact day has a minor
    # effect on sun position, but doesn't change the outcome.
    event_utc_time = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)

    # Candidate locations with their geographic coordinates and timezones
    locations = {
        "A. Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz": "America/Los_Angeles"},
        "B. Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz": "America/Chicago"},
        "C. St. John's, Newfoundland and Labrador": {"lat": 47.56, "lon": -52.71, "tz": "America/St_Johns"},
        "D. Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz": "America/Toronto"},
        "E. Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz": "Europe/London"},
    }

    # --- Analysis ---
    # A common approximation for the magnetic latitude of the auroral oval center:
    # Target Magnetic Latitude = 70° - (2 * Kp)
    target_mag_lat = 70 - (2 * kp_index)

    print(f"Analysis for a Kp={kp_index} event at {event_utc_time.strftime('%Y-%m-%d %H:%M UTC')}")
    print("\nThe final equation for the target auroral location is:")
    print(f"Target Magnetic Latitude = 70 - (2 * Kp_index)")
    print(f"Target Magnetic Latitude = 70 - (2 * {kp_index}) = {target_mag_lat:.2f}°\n")
    print("Now we check each location against this target.")
    print("-" * 110)
    print(f"{'Location':<45} {'Local Time':<12} {'Is Dark?':<12} {'Magnetic Lat.':<15} {'Closeness to Oval Center (degrees)':<15}")
    print("-" * 110)

    # Analyze each location
    results = []
    for name, data in locations.items():
        # 1. Calculate Local Time
        local_tz = pytz.timezone(data["tz"])
        local_time = event_utc_time.astimezone(local_tz)

        # 2. Check for Darkness
        # Create a simple location object for the astral library
        city_observer = type('Observer', (), {'latitude': data['lat'], 'longitude': data['lon'], 'elevation': 0})()
        s = sun(city_observer, date=event_utc_time.date(), tzinfo=local_tz)
        # We consider it "dark" if the time is after dusk and before dawn.
        is_dark = not (s["dawn"] < local_time < s["dusk"])

        # 3. Calculate Magnetic Latitude using the aacgmv2 library
        # We assume a standard aurora altitude of 110km.
        mag_lat, _, _ = aacgmv2.get_aacgm_coord(data["lat"], data["lon"], 110, event_utc_time)

        # 4. Store results for comparison
        results.append({
            "name": name,
            "local_time": local_time.strftime('%H:%M %Z'),
            "is_dark": "Yes" if is_dark else "No",
            "mag_lat": mag_lat,
            "diff_from_target": abs(mag_lat - target_mag_lat)
        })

    # --- Print Results ---
    # Sort by the closest to the target oval, which is the most important factor
    for res in sorted(results, key=lambda x: x['diff_from_target']):
        print(f"{res['name']:<45} {res['local_time']:<12} {res['is_dark']:<12} {res['mag_lat']:.2f}°{'':<10} {res['diff_from_target']:.2f}")

    print("-" * 110)
    print("\nConclusion: The best location is dark and has the smallest 'Closeness to Oval Center' value.")
    print("This indicates its magnetic latitude is nearest to the peak auroral activity for Kp=7.")

solve_aurora_location()
<<<C>>>