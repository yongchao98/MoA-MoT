def solve_aurora_location():
    """
    Analyzes which location is most likely to see an overhead aurora
    for a Kp=7 event at 06:30 UTC.
    """
    # Data for the locations including name, UTC offset, and approximate magnetic latitude.
    # In early November, North American time zones are on Standard Time.
    locations = {
        'A': {
            'name': "Portland, Oregon",
            'utc_offset': -8,
            'mag_lat': 52,
        },
        'B': {
            'name': "Madison, Wisconsin",
            'utc_offset': -6,
            'mag_lat': 54,
        },
        'C': {
            'name': "St. John's, Newfoundland and Labrador",
            'utc_offset': -3.5,
            'mag_lat': 58,
        },
        'D': {
            'name': "Alert, Nunavut",
            'utc_offset': -5,
            'mag_lat': 89,
        },
        'E': {
            'name': "Thurso, Scotland",
            'utc_offset': 0,
            'mag_lat': 61,
        }
    }

    event_utc_hr = 6
    event_utc_min = 30
    event_utc_decimal = event_utc_hr + event_utc_min / 60.0

    kp_index = 7
    # For a Kp=7 storm, overhead auroras are typically seen around 56-60 degrees magnetic latitude.
    optimal_mag_lat_min = 56
    optimal_mag_lat_max = 60

    print(f"Analysis for an overhead aurora during a Kp={kp_index} event at {event_utc_hr:02d}:{event_utc_min:02d} UTC.\n")
    print(f"The ideal magnetic latitude for an overhead aurora at Kp={kp_index} is between {optimal_mag_lat_min}° and {optimal_mag_lat_max}°.")
    print("The most intense activity often occurs in the hours around local midnight (approx. 22:00-03:00).\n")

    best_candidate = None
    max_score = -1

    # Iterate through and evaluate each location
    for key, data in locations.items():
        name = data['name']
        utc_offset = data['utc_offset']
        mag_lat = data['mag_lat']

        # 1. Calculate Local Time
        local_time_decimal = event_utc_decimal + utc_offset
        if local_time_decimal < 0:
            local_time_decimal += 24
        local_hr = int(local_time_decimal)
        local_min = int((local_time_decimal - local_hr) * 60)

        # 2. Analyze the conditions and assign a score
        print(f"--- Analysis for ({key}) {name} ---")
        
        # All calculated times are at night, so darkness is confirmed for all.
        print(f"Local Time is {local_hr:02d}:{local_min:02d}, which is dark in early November.")

        # Magnetic Latitude Analysis
        mag_lat_score = 0
        if mag_lat > optimal_mag_lat_max + 5:
            print(f"Magnetic Latitude ({mag_lat}°) is too far north (in the polar cap).")
            mag_lat_score = 1
        elif mag_lat > optimal_mag_lat_max:
            print(f"Magnetic Latitude ({mag_lat}°) is good, but slightly north of the optimal zone.")
            mag_lat_score = 4
        elif optimal_mag_lat_min <= mag_lat <= optimal_mag_lat_max:
            print(f"Magnetic Latitude ({mag_lat}°) is excellent, located within the optimal Kp={kp_index} zone.")
            mag_lat_score = 5
        else:
            print(f"Magnetic Latitude ({mag_lat}°) is too far south for an overhead view.")
            mag_lat_score = 2

        # Local Time Analysis (as a proxy for Magnetic Local Time)
        time_score = 0
        if 22 <= local_hr or local_hr < 3:
            print("Local time is in the prime window (22:00-03:00) for intense auroral activity.")
            time_score = 5
        elif 3 <= local_hr < 6:
            print("Local time is in the post-midnight/dawn sector, where aurora is still common.")
            time_score = 4
        else: 
            print("Local time is in the dawn or evening sector, less ideal for peak activity.")
            time_score = 3
        
        # Combine scores to find the best candidate
        total_score = mag_lat_score + time_score
        print(f"Resulting Likelihood Score: {total_score} out of 10")
        print("-" * 35 + "\n")

        if total_score > max_score:
            max_score = total_score
            best_candidate = f"({key}) {name}"

    print(f"Conclusion: Based on the combined analysis of magnetic latitude and local time, the location most likely to witness an overhead aurora is {best_candidate}.")

solve_aurora_location()