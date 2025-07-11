# First, you may need to install the necessary library:
# pip install aacgmv2

import datetime

def solve_aurora_location():
    """
    Analyzes locations to find the best candidate for overhead aurora viewing
    during a Kp=7 event at a specific UTC time.
    """
    # Define the event time in UTC
    event_utc = datetime.datetime(2023, 11, 5, 6, 30, 0)

    # Define the locations with their geographic coordinates [latitude, longitude]
    # and standard time zone offset from UTC.
    locations = {
        "A": {"name": "Portland, Oregon", "coords": (45.52, -122.68), "tz": -8},
        "B": {"name": "Madison, Wisconsin", "coords": (43.07, -89.40), "tz": -6},
        "C": {"name": "St. John's, Newfoundland and Labrador", "coords": (47.56, -52.71), "tz": -3.5},
        "D": {"name": "Alert, Nunavut", "coords": (82.50, -62.35), "tz": -5},
        "E": {"name": "Thurso, Scotland", "coords": (58.59, -3.52), "tz": 0}
    }

    # Ideal parameters for an overhead aurora during a Kp=7 storm
    # The main band of the aurora is typically between 55-65 degrees MLAT.
    ideal_mlat_range = (55, 65)
    # The most active part of the oval is around magnetic midnight (21:00 - 03:00 MLT).
    ideal_mlt_range = (21, 3)

    # We need the aacgmv2 library to perform the coordinate conversion.
    try:
        import aacgmv2
    except ImportError:
        print("Error: The 'aacgmv2' library is required. Please install it using 'pip install aacgmv2'")
        return

    results = []
    print("Analysis of Aurora Viewing Conditions at 06:30 UTC (Kp=7 Event)")
    print("-" * 95)
    print(f"{'Choice':<8} {'Location':<40} {'Local Time':<12} {'Mag. Lat.':<12} {'Mag. L.T.':<12} {'Score'}")
    print("-" * 95)

    for choice, data in locations.items():
        lat, lon = data["coords"]
        
        # 1. Calculate Magnetic Latitude (MLAT) and Magnetic Local Time (MLT)
        # We use an altitude of 110 km, typical for auroras.
        mlat, _, mlt = aacgmv2.get_aacgm_coord(lat, lon, 110, event_utc)
        
        # 2. Calculate local time and check for darkness
        local_time = event_utc + datetime.timedelta(hours=data["tz"])
        # Heuristic: Dark enough if between 7 PM and 7 AM local time.
        # Alert is in polar night in November, so it's always dark.
        is_dark = (local_time.hour >= 19 or local_time.hour < 7) or choice == "D"

        # 3. Score the location based on the ideal parameters
        score = 0
        if is_dark:
            # Score based on Geomagnetic Latitude
            # Being directly in the ideal range is best.
            if ideal_mlat_range[0] <= mlat <= ideal_mlat_range[1]:
                score += 2
            # Being close to the range is good.
            elif ideal_mlat_range[0] - 5 < mlat < ideal_mlat_range[1] + 5:
                score += 1
            # Being in the polar cap (too far north) is not ideal for an overhead view
            # of the main oval, which has expanded southward.
            if mlat > 75:
                score = 0
            
            # Score based on Magnetic Local Time
            # Check if MLT is in the ideal range around midnight.
            if ideal_mlt_range[0] <= mlt or mlt < ideal_mlt_range[1]:
                score += 2
        
        results.append({"choice": choice, "name": data["name"], "score": score})
        
        # Print the analysis for each location
        print(f"{choice:<8} {data['name']:<40} {local_time.strftime('%H:%M'):<12} {f'{mlat:.1f}°':<12} {f'{mlt:.1f} h':<12} {score}")

    # Find the location with the highest score
    best_location = max(results, key=lambda x: x["score"])
    
    print("-" * 95)
    print(f"\nConclusion: The location with the highest score, and thus most likely to see overhead auroras, is:")
    print(f"{best_location['choice']}. {best_location['name']}")

solve_aurora_location()