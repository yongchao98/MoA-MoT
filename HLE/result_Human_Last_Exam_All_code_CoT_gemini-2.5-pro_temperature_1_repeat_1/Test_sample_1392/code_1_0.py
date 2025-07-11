import datetime
import aacgmv2
import pytz

def find_best_aurora_location():
    """
    Calculates magnetic coordinates for several locations to determine
    the best chance of seeing overhead auroras during a Kp=7 event.
    """
    # Define the locations with their geographic coordinates (lat, lon)
    locations = {
        "A. Portland, Oregon": (45.52, -122.67),
        "B. Madison, Wisconsin": (43.07, -89.40),
        "C. St. John's, Newfoundland and Labrador": (47.56, -52.71),
        "D. Alert, Nunavut": (82.50, -62.35),
        "E. Thurso, Scotland": (58.59, -3.52)
    }

    # Define the specific time of the event in UTC
    # Year/month are representative for season (early November)
    event_time_utc = datetime.datetime(2023, 11, 5, 6, 30, tzinfo=pytz.utc)

    print(f"Analysis for aurora viewing at {event_time_utc.strftime('%Y-%m-%d %H:%M UTC')}")
    print("="*70)
    print(f"{'Location':<45} {'Mag Lat (°N)':<15} {'Mag Local Time':<15}")
    print("-"*70)

    results = []
    for name, (lat, lon) in locations.items():
        # Convert geographic coordinates to Altitude-Adjusted Corrected Geomagnetic (AACGM) coordinates
        # We use an altitude of 300 km, typical for auroral displays
        mag_lat, mag_lon, mag_lt = aacgmv2.get_aacgm_coord(lat, lon, 300, event_time_utc)
        
        # Calculate approximate standard local time for a quick check of darkness
        # Timezone offset is longitude / 15 degrees per hour
        tz_offset = datetime.timedelta(hours=lon / 15)
        local_time = event_time_utc + tz_offset
        
        print(f"{name:<45} {mag_lat:<15.2f} {mag_lt:<15.2f} (Local Time ~{local_time.strftime('%H:%M')})")
        results.append((name, mag_lat, mag_lt))

    print("="*70)
    print("\nAnalysis:")
    print("A Kp=7 event means the auroral oval is bright and expanded, typically covering 55-65° magnetic latitude.")
    print("The most intense activity occurs near magnetic midnight (approx. 22:00-02:00 MLT).")
    print("\nConclusion:")
    print("Portland and Madison are at the southern edge of the potential activity.")
    print("Alert is too far north (in the polar cap), so the main aurora would be to its south.")
    print("Thurso has a good magnetic latitude but is in the morning sector (MLT > 6), where auroras are less intense and dawn is approaching.")
    print("St. John's has an ideal magnetic latitude (58.65°) to be under the main oval and a favorable magnetic local time (2.95), placing it in the active night-side region.")

# Execute the function to print the analysis
find_best_aurora_location()
<<<C>>>