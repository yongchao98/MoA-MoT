import datetime
import pytz
# The aacgmv2 library is used for converting geographic to magnetic coordinates.
# If you don't have it, you can install it via pip:
# pip install aacgmv2
import aacgmv2

def analyze_aurora_location(utc_time, locations):
    """
    Analyzes a list of locations to determine the most likely place to see
    an overhead aurora at a specific UTC time.

    Prints the analysis for each location and a final conclusion.
    """

    # Kp=7 corresponds to an auroral oval extending down to roughly 50.4 degrees magnetic latitude.
    kp_7_boundary_mlat = 50.4

    print(f"Analysis for UTC: {utc_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Geomagnetic storm level: Kp=7 (equatorward boundary ~{kp_7_boundary_mlat}° magnetic latitude)\n")
    print("-" * 70)
    print(f"{'Location':<35} | {'Local Time':<12} | {'Mag. Lat.':<12} | {'Mag. L.T.':<12}")
    print("-" * 70)

    results = []

    for name, data in locations.items():
        # Calculate local time
        local_tz = pytz.timezone(data['tz'])
        local_time = utc_time.astimezone(local_tz)

        # Calculate Altitude-Adjusted Corrected Geo-Magnetic (AACGM) coordinates
        # We use an altitude of 300 km, typical for auroral phenomena.
        mlat, mlon, mlt = aacgmv2.get_aacgm_coord(data['lat'], data['lon'], 300, utc_time)

        results.append({
            "name": name,
            "local_time_str": local_time.strftime('%H:%M'),
            "mlat": mlat,
            "mlt": mlt,
        })
        
        # Print results for each location
        print(f"{name:<35} | {local_time.strftime('%H:%M %Z'):<12} | {mlat:^12.1f} | {mlt:^12.2f}")

    print("-" * 70)

    # --- Analysis for the conclusion ---
    best_candidate = None
    max_score = -1

    for res in results:
        # Scoring criteria:
        # 1. Must be dark (roughly 19:00 - 05:00 local time).
        # 2. Must be within the auroral zone (mlat > Kp 7 boundary).
        # 3. Higher magnetic latitude is better (but not in polar cap, > ~80).
        # 4. MLT should be close to midnight (21:00-03:00).

        # Check darkness
        hour = int(res['local_time_str'].split(':')[0])
        is_dark = 19 <= hour <= 23 or 0 <= hour <= 5

        # Check magnetic latitude
        in_zone = res['mlat'] > kp_7_boundary_mlat
        in_polar_cap = res['mlat'] > 80

        # Check Magnetic Local Time
        in_prime_mlt = 21 <= res['mlt'] <= 24 or 0 <= res['mlt'] <= 3

        # Simple scoring
        score = 0
        if is_dark and in_zone and not in_polar_cap:
            # Score bonus for being deep in the oval
            score += res['mlat']
            # Score bonus for being close to magnetic midnight
            if in_prime_mlt:
                score += 10 # Add a significant bonus for prime time

        if score > max_score:
            max_score = score
            best_candidate = res['name']
            
    print("\nConclusion:")
    print("To be the *most likely* location for an *overhead* aurora, a site needs to:")
    print("1. Be in darkness.")
    print("2. Be located deep inside the auroral oval (high magnetic latitude).")
    print("3. Be near magnetic midnight (peak activity time).")
    print("\n- Alert is in the polar cap, too far north for the main oval.")
    print("- Thurso is in its local morning, so it is not in darkness and is outside the prime night-side viewing sector.")
    print("- Portland and Madison are at lower magnetic latitudes, near the edge of the oval for a Kp=7 storm. Auroras might be visible, but likely on the horizon, not overhead.")
    print("- St. John's has a high magnetic latitude (~57°), placing it firmly inside the main auroral oval. Even though its magnetic local time is past midnight, its favorable latitude makes it the most probable location to be directly under the aurora.")

    print(f"\nTherefore, the most likely location is: {best_candidate}")


if __name__ == '__main__':
    # Define the time of the event
    # Using a representative date in early November. Time is 06:30 UTC.
    utc_time_event = datetime.datetime(2023, 11, 5, 6, 30, 0, tzinfo=pytz.utc)

    # Define the locations with their geographic coordinates and time zones
    locations_data = {
        "A. Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz": "America/Los_Angeles"}, # UTC-8
        "B. Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz": "America/Chicago"}, # UTC-6
        "C. St. John's, Newfoundland": {"lat": 47.56, "lon": -52.71, "tz": "America/St_Johns"}, # UTC-3.5
        "D. Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz": "America/Toronto"}, # UTC-5
        "E. Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz": "Europe/London"}, # UTC+0
    }

    analyze_aurora_location(utc_time_event, locations_data)
<<<C>>>