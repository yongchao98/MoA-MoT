import datetime
import pytz
from timezonefinder import TimezoneFinder
try:
    import aacgmv2
except ImportError:
    print("Required library not found. Please install it by running: pip install aacgmv2")
    exit()

def solve_aurora_location():
    """
    Analyzes locations to find the most likely spot for an overhead aurora
    during a Kp=7 event at a specific UTC time.
    """
    # Define the locations with their geographic coordinates (latitude, longitude)
    locations = {
        "A. Portland, Oregon": (45.52, -122.68),
        "B. Madison, Wisconsin": (43.07, -89.40),
        "C. St. John's, Newfoundland and Labrador": (47.56, -52.71),
        "D. Alert, Nunavut": (82.50, -62.35),
        "E. Thurso, Scotland": (58.59, -3.52)
    }

    # Define the event time and Kp index
    # The specific date doesn't matter, only the time of day for lighting conditions.
    event_utc_time = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)
    kp = 7

    # --- Step 1: Estimate the auroral oval's central magnetic latitude ---
    # A common approximation for the center of the auroral oval is based on the Kp index.
    # Formula: Center Magnetic Latitude = 66.5 - 1.5 * Kp
    print("--- Auroral Oval Position Analysis ---")
    print("The aurora is most likely to be overhead near the center of the auroral oval.")
    print("We can estimate the magnetic latitude of the oval's center with the formula: 66.5 - 1.5 * Kp")
    
    # Outputting each number in the final equation as requested
    val1, val2, val3 = 66.5, 1.5, kp
    auroral_center_lat = val1 - val2 * val3
    print(f"For Kp={kp}, the calculation is: {val1} - {val2} * {val3} = {auroral_center_lat:.1f}°\n")

    # --- Step 2: Analyze each location ---
    print("--- Location Analysis ---")
    print(f"Analyzing each location for conditions at {event_utc_time.strftime('%H:%M')} UTC:")
    
    tf = TimezoneFinder()
    results = []

    for name, (lat, lon) in locations.items():
        # Calculate local time to check for darkness
        tz_name = tf.timezone_at(lng=lon, lat=lat)
        local_tz = pytz.timezone(tz_name)
        local_time = event_utc_time.astimezone(local_tz)
        
        # A simple check for darkness (e.g., between 8 PM and 5 AM)
        is_dark = (local_time.hour >= 20 or local_time.hour < 5)

        # Calculate magnetic latitude using the aacgmv2 library
        # We use an altitude of 300km, a typical altitude for auroras.
        mlat, _, _ = aacgmv2.get_aacgm_coord(lat, lon, 300, event_utc_time)

        # Calculate how close the location's magnetic latitude is to the oval's center
        diff_from_oval = abs(mlat - auroral_center_lat)

        results.append({
            "name": name,
            "local_time": local_time.strftime('%H:%M %Z'),
            "is_dark": is_dark,
            "mag_lat": mlat,
            "diff": diff_from_oval
        })

    # --- Step 3: Print results and find the best candidate ---
    print(f"\n{'Location':<40} {'Local Time':<12} {'Dark?':<7} {'Magnetic Lat.':<15} {'Proximity to Oval Center':<25}")
    print("-" * 110)

    # Sort by proximity to the oval, but we will make the final decision based on darkness as well
    sorted_results = sorted(results, key=lambda x: x['diff'])
    
    best_candidate = None
    min_diff = float('inf')

    for res in sorted_results:
        print(f"{res['name']:<40} {res['local_time']:<12} {str(res['is_dark']):<7} {res['mag_lat']:.2f}°{'':<8} {res['diff']:.2f}°")
        # The best candidate must be dark
        if res['is_dark'] and res['diff'] < min_diff:
            min_diff = res['diff']
            best_candidate = res['name']

    print("\n--- Conclusion ---")
    print("An overhead aurora requires the location to be both in darkness and have a magnetic latitude")
    print(f"very close to the estimated auroral oval center of {auroral_center_lat:.1f}°.")
    print(f"\nThe location that best fits these criteria is: {best_candidate}")

solve_aurora_location()
<<<C>>>