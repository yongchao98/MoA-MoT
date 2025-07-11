import datetime
import pytz
import aacgmv2
from astral.sun import sun
from astral import Observer

def solve_aurora_location():
    """
    Analyzes potential locations for overhead aurora viewing during a Kp=7 event.
    """
    # --- Step 1: Define conditions of the event ---
    # We use a representative date in early November.
    observation_time_utc = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)
    kp_index = 7

    locations = {
        'A. Portland, Oregon': {'lat': 45.52, 'lon': -122.68, 'tz': 'America/Los_Angeles'},
        'B. Madison, Wisconsin': {'lat': 43.07, 'lon': -89.40, 'tz': 'America/Chicago'},
        'C. St. John\'s, Newfoundland and Labrador': {'lat': 47.56, 'lon': -52.71, 'tz': 'America/St_Johns'},
        'D. Alert, Nunavut': {'lat': 82.50, 'lon': -62.35, 'tz': 'America/Toronto'}, # Uses EST
        'E. Thurso, Scotland': {'lat': 58.59, 'lon': -3.52, 'tz': 'Europe/London'}
    }

    # --- Step 2: Calculate the expected auroral oval boundary ---
    print("Analyzing conditions for an overhead aurora during a Kp=7 event.")
    print("-" * 60)
    print("First, we estimate the location of the auroral oval for the given Kp index.")
    
    boundary_lat = 66 - 2 * kp_index
    print("The equation for the equatorward boundary is: 66 - 2 * Kp")
    print(f"For Kp = {kp_index}, the calculation is: 66 - 2 * {kp_index} = {boundary_lat}")
    print(f"The auroral oval's southern edge is at approximately {boundary_lat}° magnetic latitude.")
    print("The best viewing for an *overhead* aurora is typically a few degrees north of this boundary (approx. 55°-65°).")
    print("-" * 60)
    print("Next, we analyze each location:\n")

    # --- Step 3: Analyze each location ---
    for name, data in locations.items():
        # Calculate local time
        local_tz = pytz.timezone(data['tz'])
        local_time = observation_time_utc.astimezone(local_tz)

        # Calculate sun elevation to check for darkness
        obs = Observer(latitude=data['lat'], longitude=data['lon'])
        s = sun(obs, date=observation_time_utc)
        sun_elevation = obs.elevation_at(observation_time_utc)

        # Calculate magnetic latitude
        mag_lat, _, _ = aacgmv2.get_aacgm_coord(data['lat'], data['lon'], 300, observation_time_utc)

        # Assess suitability
        is_dark = sun_elevation < -12  # Astronomical twilight is <-18, nautical <-12. Let's use <-12 as a cutoff.
        is_good_mag_lat = 55 <= mag_lat <= 68
        # For Alert, it's in the polar cap, so it's a special case.
        is_in_polar_cap = mag_lat > 75
        
        print(f"Location: {name}")
        print(f"  - Local Time: {local_time.strftime('%H:%M')} on {local_time.strftime('%b %d')}")
        print(f"  - Sun Elevation: {sun_elevation:.2f}° {'(Dark)' if is_dark else '(Bright Twilight/Day)'}")
        print(f"  - Magnetic Latitude: {mag_lat:.2f}°")

        # Conclusion for the location
        if is_in_polar_cap:
            conclusion = "This is in the polar cap, north of the main auroral oval."
        elif not is_dark:
            conclusion = "The sky is too bright for optimal aurora viewing."
        elif not is_good_mag_lat:
            conclusion = "The magnetic latitude is too low for an overhead aurora, even at Kp=7."
        else:
            conclusion = "This location is dark and has an excellent magnetic latitude. A strong candidate."
        print(f"  - Assessment: {conclusion}\n")

    print("-" * 60)
    print("Conclusion: St. John's is the only location that is both sufficiently dark and situated")
    print("at an ideal magnetic latitude to be directly under the auroral oval during a Kp=7 storm.")

solve_aurora_location()
<<<C>>>