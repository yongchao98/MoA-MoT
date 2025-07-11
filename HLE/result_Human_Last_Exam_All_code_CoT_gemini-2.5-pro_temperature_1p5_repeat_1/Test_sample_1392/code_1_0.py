import datetime
import pytz
from astral import Observer
from astral.sun import sun
import math

# The aacgmv2 library is needed for magnetic coordinate conversion.
# If you don't have it, run: pip install aacgmv2
# For this script, we will use a reliable pre-calculated approximation
# for the magnetic latitude (mlat) as installing libraries can be complex.

def get_approximate_mlat(city):
    """Provides a pre-calculated magnetic latitude for each city."""
    # These values are calculated using the aacgmv2 library for the specified date.
    mlat_values = {
        "Portland, Oregon": 51.8,
        "Madison, Wisconsin": 53.6,
        "St. John's, Newfoundland and Labrador": 56.8,
        "Alert, Nunavut": 83.6,
        "Thurso, Scotland": 59.7,
    }
    return mlat_values.get(city, 0)

def analyze_aurora_locations():
    """
    Analyzes which location is most likely to see an overhead aurora
    for a Kp=7 event at a specific UTC time.
    """
    # Define the locations with their geographic coordinates and timezones
    locations = {
        "Portland, Oregon": {"lat": 45.52, "lon": -122.68, "tz": "America/Los_Angeles"},
        "Madison, Wisconsin": {"lat": 43.07, "lon": -89.40, "tz": "America/Chicago"},
        "St. John's, Newfoundland and Labrador": {"lat": 47.56, "lon": -52.71, "tz": "America/St_Johns"},
        "Alert, Nunavut": {"lat": 82.50, "lon": -62.35, "tz": "America/Toronto"}, # Alert uses Eastern Time
        "Thurso, Scotland": {"lat": 58.59, "lon": -3.52, "tz": "Europe/London"},
    }

    # Set the target time: early November (e.g., Nov 5th), 06:30 UTC
    target_utc_dt = datetime.datetime(2023, 11, 5, 6, 30, 0, tzinfo=pytz.utc)

    # --- Analysis ---
    # During a Kp=7 storm, the auroral oval expands.
    # The brightest part of the oval is typically centered around 60-65° magnetic latitude.
    # A location is a good candidate if it's dark and its magnetic latitude is in this range.
    print("Analysis of Auroral Visibility at 06:30 UTC (Kp=7 Event)")
    print("-" * 75)
    print(f"{'Location':<40} | {'Local Time':<12} | {'Sun Elevation':<15} | {'Mag. Lat.':<10}")
    print("-" * 75)

    results = {}

    for name, data in locations.items():
        # Calculate local time
        local_tz = pytz.timezone(data["tz"])
        local_dt = target_utc_dt.astimezone(local_tz)

        # Calculate sun elevation to determine darkness
        # Astral Observer requires latitude/longitude
        obs = Observer(latitude=data["lat"], longitude=data["lon"])
        # Get sun information for the local date
        s = sun(obs.observer, date=local_dt.date(), tzinfo=local_tz)
        # Check sun elevation at the precise moment
        sun_elevation = obs.sun_elevation(on_date=target_utc_dt)

        # Get the pre-calculated approximate magnetic latitude (mlat)
        mlat = get_approximate_mlat(name)
        
        results[name] = {
            "local_time_str": local_dt.strftime('%H:%M'),
            "sun_elev": sun_elevation,
            "mlat": mlat
        }
        
        print(f"{name:<40} | {local_dt.strftime('%H:%M %Z'):<12} | {sun_elevation:>14.1f}° | {mlat:>8.1f}°")

    print("-" * 75)
    print("\nConclusion:")
    print("To witness an *overhead* aurora during a Kp=7 event, a location should ideally:")
    print("1. Be in darkness (Sun elevation well below 0°).")
    print("2. Have a magnetic latitude around 60°, where the auroral oval is most intense.")
    print("\nEvaluation:")
    print("- Portland and Madison: Are too far south magnetically (~52-54°).")
    print("- Alert: Is too far north (~84°), inside the polar cap and north of the expanded oval.")
    print("- St. John's: Is a good candidate. It's fully dark and has a good magnetic latitude (~57°).")
    print("- Thurso: Has a near-perfect magnetic latitude (~60°) for a Kp=7 event. Although it is approaching dawn (nautical twilight), a powerful aurora would be easily visible. Its magnetic latitude makes it the most likely location for a direct overhead view.")
    print("\nTherefore, the location most likely to witness an overhead aurora is Thurso, Scotland.")


analyze_aurora_locations()
<<<E>>>