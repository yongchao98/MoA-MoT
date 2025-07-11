import math

def calculate_analysis():
    """
    Analyzes which location is most likely to see overhead auroras for a given
    geomagnetic storm level (Kp-index) and time (UTC).
    """
    # --- Input Parameters ---
    # Event time in UTC (06:30)
    utc_time_hours = 6.5
    # Kp-index of the geomagnetic storm
    kp_index = 7

    # Location of the North Geomagnetic Pole (approximate value for calculation)
    # Note: The pole moves. This is a standard approximation.
    lat_gp, lon_gp = 80.7, -72.7

    # City data: Name -> (Geographic Latitude, Geographic Longitude)
    locations = {
        "Portland, Oregon": (45.5, -122.7),
        "Madison, Wisconsin": (43.1, -89.4),
        "St. John's, Newfoundland and Labrador": (47.6, -52.7),
        "Alert, Nunavut": (82.5, -62.3),
        "Thurso, Scotland": (58.6, -3.5)
    }

    print("--- Auroral Sighting Analysis ---")
    print(f"Analyzing potential for overhead aurora for a Kp={kp_index} event at {int(utc_time_hours):02}:{int((utc_time_hours % 1) * 60):02} UTC.\n")

    # Ideal auroral zone for a given Kp index. This is an empirical formula.
    # The oval's center is roughly at this magnetic latitude.
    ideal_mag_lat_eq = 66 - kp_index * 2
    # The oval has a width, so we define a "good" range.
    # For Kp=7, the equatorward edge is around 56-58 degrees. Let's use a prime viewing zone of 56-66.
    auroral_zone_low = 56
    auroral_zone_high = 66

    print(f"For Kp={kp_index}, the ideal magnetic latitude for overhead aurora is ~{auroral_zone_low}° to {auroral_zone_high}°.\n")

    best_candidate = None
    highest_score = -1

    for name, (lat, lon) in locations.items():
        print(f"--- Evaluating: {name} ---")

        # 1. Calculate Local Time
        local_time_hours = (utc_time_hours + lon / 15.0) % 24
        is_dark = local_time_hours < 5.0 or local_time_hours > 20.0
        time_str = f"{int(local_time_hours):02}:{int((local_time_hours % 1) * 60):02}"
        
        # We need to handle times on the previous day (e.g., Portland)
        if utc_time_hours + lon / 15.0 < 0:
            local_time_hours += 24
            time_str = f"{int(local_time_hours):02}:{int((local_time_hours % 1) * 60):02} (Previous Day)"


        print(f"1. Local Time Check:")
        print(f"   - Local time at {lon}° longitude is {time_str}.")
        print(f"   - Is it dark? {'Yes' if is_dark else 'No (Daylight/Twilight)'}")

        # 2. Calculate Magnetic Latitude
        lat_rad, lon_rad = math.radians(lat), math.radians(lon)
        lat_gp_rad, lon_gp_rad = math.radians(lat_gp), math.radians(lon_gp)
        
        # Using the spherical law of cosines to find the angular distance (p) to the geomagnetic pole
        cos_p = (math.sin(lat_rad) * math.sin(lat_gp_rad) +
                 math.cos(lat_rad) * math.cos(lat_gp_rad) * math.cos(lon_rad - lon_gp_rad))
        
        # Ensure the value is within the domain of acos to avoid math errors
        cos_p = max(-1.0, min(1.0, cos_p))
        
        p_rad = math.acos(cos_p)
        mag_colatitude = math.degrees(p_rad)
        mag_latitude = 90.0 - mag_colatitude
        
        # Check if it is within the prime auroral zone
        in_zone = auroral_zone_low <= mag_latitude <= auroral_zone_high

        print(f"2. Magnetic Latitude Check:")
        print(f"   - Calculation: Mag Lat = 90 - arccos(sin({lat:.1f})*sin({lat_gp:.1f}) + cos({lat:.1f})*cos({lat_gp:.1f})*cos({lon:.1f} - ({lon_gp:.1f})))")
        print(f"   - Resulting Magnetic Latitude: {mag_latitude:.1f}°")
        print(f"   - Is this in the Kp={kp_index} auroral zone ({auroral_zone_low}°-{auroral_zone_high}°)? {'Yes' if in_zone else 'No'}")

        # 3. Conclusion for this location
        score = 0
        if is_dark:
            score += 1
        if in_zone:
            score += 2 # Give more weight to being in the zone
        
        if score > highest_score:
            highest_score = score
            best_candidate = name
        
        print(f"   - Verdict: {('A promising candidate.' if score == 3 else 'Unlikely to see overhead aurora.')}\n")
    
    print("--- FINAL CONCLUSION ---")
    if best_candidate:
        print(f"The location most likely to witness overhead auroras is {best_candidate}.")
        print("This is because it is dark at the specified time and its magnetic latitude falls squarely within the expanded auroral oval for a Kp=7 storm.")
    else:
        print("No suitable candidate found based on the criteria.")

calculate_analysis()
<<<C>>>