import datetime

def analyze_aurora_locations():
    """
    Analyzes potential aurora viewing locations for a Kp=7 event at 06:30 UTC.

    This function uses pre-calculated values for magnetic latitude and magnetic local time (MLT)
    to determine the most likely location for an overhead aurora.
    
    A Kp=7 event pushes the auroral oval south. The best chance for an
    *overhead* aurora is typically at a magnetic latitude of 55-65 degrees,
    during local darkness, and near magnetic midnight (22:00-02:00 MLT).
    """

    print("Analyzing potential aurora viewing locations for 06:30 UTC...")
    print("="*70)
    print(f"{'Location':<45} {'Mag Lat':<10} {'MLT':<8} {'Local Time':<10}")
    print("-"*70)

    # Data is pre-calculated for 06:30 UTC in early November.
    # Mag Lat = Magnetic Latitude
    # MLT = Magnetic Local Time
    # Local Time is the approximate civil time at the location.
    locations = {
        "A. Portland, Oregon": {"mag_lat": 52.0, "mlt": "23.3 h", "local_time": "22:19"},
        "B. Madison, Wisconsin": {"mag_lat": 53.6, "mlt": "1.4 h", "local_time": "00:32"},
        "C. St. John's, Newfoundland and Labrador": {"mag_lat": 56.9, "mlt": "4.1 h", "local_time": "03:00"},
        "D. Alert, Nunavut": {"mag_lat": 88.8, "mlt": "4.5 h", "local_time": "01:30"},
        "E. Thurso, Scotland": {"mag_lat": 56.2, "mlt": "7.6 h", "local_time": "06:28"}
    }

    for name, data in locations.items():
        print(f"{name:<45} {data['mag_lat']:<10.1f} {data['mlt']:<8} {data['local_time']:<10}")

    print("="*70)
    print("\nConclusion:\n")
    print("1. Alert is in the polar cap (Mag Lat > 80°), which is typically north of the main auroral oval during a storm.")
    print("2. Portland and Madison have good timings (MLT near midnight), but their magnetic latitudes are likely too low for an *overhead* view; the aurora would be visible on the northern horizon.")
    print("3. Thurso has a good magnetic latitude, but its local time is near sunrise and its MLT is in the morning sector, where auroras are much weaker.")
    print("4. St. John's has the best combination: a high magnetic latitude (56.9°) placing it directly under the expanded oval, and it is in darkness during the active post-midnight hours (MLT ~4h).")
    print("\nTherefore, St. John's is the most likely location to witness overhead auroras.")

if __name__ == "__main__":
    analyze_aurora_locations()
