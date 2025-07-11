import datetime

def calculate_aurora_chances():
    """
    Analyzes which location is most likely to see overhead auroras
    during a Kp=7 event at 06:30 UTC in early November.
    """
    event_utc_time_str = "06:30"
    utc_hour, utc_minute = map(int, event_utc_time_str.split(':'))

    locations = {
        "A. Portland, Oregon": {
            "utc_offset": -8.0,
            "mag_lat": 52,
            "analysis": "Dark sky. On the southern edge of the Kp=7 auroral oval. Overhead aurora is less likely."
        },
        "B. Madison, Wisconsin": {
            "utc_offset": -6.0,
            "mag_lat": 53,
            "analysis": "Dark sky. On the southern edge of the Kp=7 auroral oval. Overhead aurora is less likely."
        },
        "C. St. John's, Newfoundland and Labrador": {
            "utc_offset": -3.5,
            "mag_lat": 56,
            "analysis": "Deep night (very dark). Excellent position under the expanded Kp-7 auroral oval. Prime candidate."
        },
        "D. Alert, Nunavut": {
            "utc_offset": -5.0,
            "mag_lat": 86,
            "analysis": "Dark sky (polar night). Too far north (inside the polar cap). Aurora would be on the southern horizon, not overhead."
        },
        "E. Thurso, Scotland": {
            "utc_offset": 0.0,
            "mag_lat": 61,
            "analysis": "Time is at sunrise (dawn). Brightening sky would wash out the aurora despite a good magnetic latitude."
        }
    }

    print(f"Analyzing aurora visibility for a Kp=7 event at {event_utc_time_str} UTC.\n")
    print("-" * 70)
    print(f"{'Location':<45} {'Local Time':<15} {'Analysis'}")
    print("-" * 70)

    for name, data in locations.items():
        offset = datetime.timedelta(hours=data["utc_offset"])
        local_time = datetime.datetime(2023, 11, 5, utc_hour, utc_minute) + offset
        local_time_str = local_time.strftime('%H:%M')
        
        # Handle cases where the day changes
        if offset.total_seconds() < 0:
             local_time_str += " (prev day)"
        
        print(f"{name:<45} {local_time_str:<15} {data['analysis']}")

    print("\n" + "="*70)
    print("Conclusion: St. John's offers the best combination of a dark sky (middle of the night) and an ideal magnetic latitude to be under the expanded auroral oval.")
    print("="*70)


calculate_aurora_chances()