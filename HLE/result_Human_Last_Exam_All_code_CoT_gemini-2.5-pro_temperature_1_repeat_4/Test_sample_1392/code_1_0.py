import datetime
try:
    import aacgmv2
except ImportError:
    print("Please install the aacgmv2 library by running: pip install aacgmv2")
    exit()

def solve_aurora_location():
    """
    Calculates the best location to see an overhead aurora for a Kp=7 event
    at a specific time by computing the magnetic latitude and magnetic local time for each option.
    """
    # Define the time of the event
    # An arbitrary date in early November is chosen as the specific day has a minor effect on the result.
    event_time_utc = datetime.datetime(2023, 11, 5, 6, 30, 0)

    # Define the locations with their geographic coordinates (Latitude, Longitude)
    locations = {
        'A': {"name": "Portland, Oregon", "coords": (45.52, -122.67)},
        'B': {"name": "Madison, Wisconsin", "coords": (43.07, -89.40)},
        'C': {"name": "St. John's, Newfoundland and Labrador", "coords": (47.56, -52.71)},
        'D': {"name": "Alert, Nunavut", "coords": (82.50, -62.35)},
        'E': {"name": "Thurso, Scotland", "coords": (58.59, -3.52)}
    }

    # Kp=7 event characteristics:
    # The auroral oval expands. Overhead views are most likely in the 58° to 65° magnetic latitude band.
    # Prime activity occurs around magnetic midnight (21:00-03:00 MLT).
    kp_index = 7
    prime_mlat_min = 58
    prime_mlat_max = 65
    prime_mlt_min = 21.0
    prime_mlt_max = 3.0 # Wraps around midnight

    print(f"Analyzing potential aurora viewing locations for a Kp={kp_index} event at {event_time_utc} UTC.\n")
    print(f"Ideal conditions: Magnetic Latitude between {prime_mlat_min}° and {prime_mlat_max}°, and Magnetic Local Time near midnight.\n")

    best_candidate = None
    best_score = -1

    for key, loc in locations.items():
        lat, lon = loc["coords"]
        # Altitude for auroral emissions is typically 100-400 km. We use an average of 300 km.
        mlat, mlon, mlt = aacgmv2.get_aacgm_coord(lat, lon, 300, event_time_utc)

        print(f"--- Option {key}: {loc['name']} ---")
        print(f"Magnetic Latitude: {mlat:.2f}°")
        print(f"Magnetic Local Time: {mlt:.2f} hrs")

        # Scoring logic to find the best candidate
        # Score is higher for being in the prime MLAT band.
        # Bonus points for being in the prime MLT window.
        score = 0
        if prime_mlat_min <= mlat <= prime_mlat_max:
            score += 10  # Prime candidate based on latitude
            print("Analysis: This location is in the prime magnetic latitude band for an overhead aurora.")
        elif mlat > prime_mlat_max:
            print("Analysis: This location is too far north (in the polar cap), less likely to see the main auroral oval overhead.")
        elif mlat > 50: # Plausible viewing, but maybe not overhead
            print("Analysis: This location is at a plausible latitude, but may be on the southern edge of the oval.")
        else:
            print("Analysis: This location is too far south.")

        # Check for Magnetic Local Time
        # Handle wraparound time (e.g., 21:00 to 03:00)
        in_prime_time = False
        if prime_mlt_min > prime_mlt_max: # e.g. 21:00 to 03:00
            if mlt >= prime_mlt_min or mlt <= prime_mlt_max:
                in_prime_time = True
        else: # e.g. 09:00 to 17:00
            if prime_mlt_min <= mlt <= prime_mlt_max:
                in_prime_time = True
        
        if in_prime_time:
            score += 5 # Bonus for prime timing
            print("Analysis: The timing is excellent (near magnetic midnight).")
        else:
             print("Analysis: The timing is outside the peak activity window (magnetic midnight).")
        
        if score > best_score:
            best_score = score
            best_candidate = key
        
        print("\n")


    print("="*40)
    print("Conclusion:")
    print(f"Based on the analysis, the location with the highest likelihood of seeing an overhead aurora is Option {best_candidate}: {locations[best_candidate]['name']}.")
    print("This is because its magnetic latitude falls squarely within the expected range for a Kp=7 auroral oval, and its magnetic local time is well within the window of peak auroral activity.")

solve_aurora_location()
<<<C>>>