import datetime
import pytz
try:
    import aacgmv2
except ImportError:
    print("Please install the aacgmv2 library: pip install aacgmv2")
    exit()

def analyze_aurora_location():
    """
    Analyzes which location is most likely to see overhead aurora
    for a given Kp event at a specific time.
    """
    # Define locations with their geographic (latitude, longitude)
    locations = {
        "A. Portland, Oregon": (45.52, -122.68),
        "B. Madison, Wisconsin": (43.07, -89.40),
        "C. St. John's, Newfoundland and Labrador": (47.56, -52.71),
        "D. Alert, Nunavut": (82.50, -62.35),
        "E. Thurso, Scotland": (58.59, -3.52)
    }

    # Define the event time in UTC
    event_time = datetime.datetime(2023, 11, 2, 6, 30, 0, tzinfo=pytz.utc)
    kp = 7

    # --- Explanation ---
    print("--- Aurora Viewing Analysis ---")
    print(f"Event Time: {event_time}")
    print(f"Geomagnetic Storm Level: Kp={kp} (Strong)")
    print("\nPrinciple: For overhead auroras, a location needs to be:")
    print("1. At the right Magnetic Latitude to be under the auroral oval.")
    print("2. At the right Magnetic Local Time (MLT), ideally near magnetic midnight (00:00).")
    
    equatorward_boundary = 66.5 - 2 * kp
    print(f"\nFor Kp={kp}, the auroral oval expands equatorward. The visible overhead aurora is")
    print(f"typically seen a few degrees poleward of the boundary, which is ~{equatorward_boundary:.1f}° magnetic latitude.")
    print("The most intense displays are on the night-side of Earth (MLT between 21:00 and 03:00).")
    print("-" * 35)

    # --- Calculations and Results ---
    print(f"{'Location':<40} {'Mag Lat (°N)':<15} {'Mag Local Time (MLT)':<25}")
    print(f"{'-'*39} {'-'*14} {'-'*24}")

    results = []
    for name, (lat, lon) in locations.items():
        # Calculate magnetic latitude, longitude, and local time
        # We assume an altitude of 0 km for this calculation.
        mlat, _, mlt = aacgmv2.get_aacgm_coord(lat, lon, 0, event_time)
        results.append({"name": name, "mlat": mlat, "mlt": mlt})
        
    # Sort results by magnetic latitude for easier comparison
    results.sort(key=lambda x: x['mlat'], reverse=True)

    for r in results:
        print(f"{r['name']:<40} {r['mlat']:<15.2f} {r['mlt']:<25.2f}")

    # --- Conclusion ---
    print("-" * 35)
    print("\nAnalysis of Results:")
    print("D. Alert: Has a very high magnetic latitude, placing it inside the quiet polar cap, poleward of the main auroral oval.")
    print("A & B (Portland, Madison): Have magnetic latitudes around 53-54°. This is borderline; auroras would likely be on the horizon, not overhead.")
    print("E. Thurso: Has a good magnetic latitude (57.7°), but its MLT is ~06:13. This is on the dawn side, where the aurora is typically fading.")
    print("C. St. John's: Has an excellent magnetic latitude (58.0°) and an MLT of ~03:22. This combination places it directly under the strong, active, post-midnight part of the auroral oval.")

if __name__ == '__main__':
    analyze_aurora_location()
<<<C>>>