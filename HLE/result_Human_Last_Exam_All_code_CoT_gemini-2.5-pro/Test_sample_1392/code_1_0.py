import datetime

def solve_aurora_location():
    """
    Analyzes locations to determine the most likely spot to see an overhead aurora
    during a Kp=7 event at 06:30 UTC in early November.
    """

    # Data for each location:
    # 'name': Location name
    # 'lat': Geographic Latitude
    # 'lon': Geographic Longitude
    # 'utc_offset': UTC offset in hours
    # 'mag_lat': Approximate Magnetic Latitude
    locations = {
        'A': {'name': 'Portland, Oregon', 'utc_offset': -8, 'mag_lat': 52},
        'B': {'name': 'Madison, Wisconsin', 'utc_offset': -6, 'mag_lat': 53},
        'C': {'name': "St. John's, Newfoundland and Labrador", 'utc_offset': -3.5, 'mag_lat': 58},
        'D': {'name': 'Alert, Nunavut', 'utc_offset': -5, 'mag_lat': 86},
        'E': {'name': 'Thurso, Scotland', 'utc_offset': 0, 'mag_lat': 61},
    }

    # Event details
    event_utc_hour = 6
    event_utc_minute = 30
    kp_index = 7
    # For Kp=7, the auroral oval is intense and expands equatorward.
    # The ideal magnetic latitude for an *overhead* view is roughly 55-60 degrees.
    # Locations with magnetic latitudes from ~50 to ~65 degrees are plausible.
    prime_mag_lat_min = 55
    prime_mag_lat_max = 60

    print(f"Analyzing likelihood of overhead aurora for a Kp={kp_index} event at {event_utc_hour:02d}:{event_utc_minute:02d} UTC.\n")
    print("The ideal magnetic latitude for an overhead view is considered to be between 55° and 60°.")
    print("The local time must also be dark for viewing.\n")

    best_candidate = None
    highest_score = -1

    for key, loc in locations.items():
        # Calculate local time
        local_time_h = event_utc_hour + loc['utc_offset']
        local_time_m = event_utc_minute
        if local_time_h < 0:
            local_time_h += 24
        if local_time_h >= 24:
            local_time_h -= 24
        
        # Assess suitability
        # Score based on magnetic latitude and darkness
        score = 0
        
        # 1. Magnetic Latitude Check
        mag_lat = loc['mag_lat']
        mag_lat_comment = ""
        if mag_lat >= prime_mag_lat_min and mag_lat <= prime_mag_lat_max:
            mag_lat_comment = f"Excellent magnetic latitude ({mag_lat}°), right in the prime zone."
            score += 2
        elif mag_lat > 50 and mag_lat < 65:
            mag_lat_comment = f"Good magnetic latitude ({mag_lat}°), on the edge of the prime zone."
            score += 1
        elif mag_lat > 65:
            mag_lat_comment = f"Poor magnetic latitude ({mag_lat}°), likely too far north (in the polar cap)."
            score -= 1
        else:
             mag_lat_comment = f"Poor magnetic latitude ({mag_lat}°), likely too far south."
             score -= 1


        # 2. Darkness Check (Early November)
        # 20:00 (8 PM) to 05:00 (5 AM) is a great viewing window. 05:00-07:00 is pre-dawn, still possible.
        local_time_str = f"{int(local_time_h):02d}:{local_time_m:02d}"
        darkness_comment = ""
        if 20 <= local_time_h < 24 or 0 <= local_time_h < 5:
            darkness_comment = f"Excellent viewing time ({local_time_str} local), sky is very dark."
            score += 2
        elif 5 <= local_time_h < 7 or 18 <= local_time_h < 20:
            darkness_comment = f"Fair viewing time ({local_time_str} local), during twilight hours."
            score += 1
        else:
            darkness_comment = f"Poor viewing time ({local_time_str} local), likely daylight."
            score -= 1
        
        print(f"--- Analysis for ({key}) {loc['name']} ---")
        print(f"Result: {mag_lat_comment}")
        print(f"Result: {darkness_comment}")
        print("-" * (len(loc['name']) + 21))
        print("")

        if score > highest_score:
            highest_score = score
            best_candidate = key

    print(f"\nConclusion: Based on the analysis, the location with the best combination of magnetic latitude and local darkness is ({best_candidate}).")


solve_aurora_location()
<<<C>>>