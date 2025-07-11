import datetime

def analyze_aurora_locations():
    """
    Analyzes which location is most likely to see overhead auroras for a given event.
    """
    event_utc_hour = 6
    event_utc_minute = 30
    kp_index = 7

    # Data for each location: {name: (geo_lat, geo_lon, utc_offset, mag_lat)}
    # Magnetic latitudes are approximate but sufficient for this comparison.
    locations = {
        "A. Portland, Oregon": {"lat": 45.5, "lon": -122.7, "utc_offset": -8.0, "mag_lat": 52.0},
        "B. Madison, Wisconsin": {"lat": 43.1, "lon": -89.4, "utc_offset": -6.0, "mag_lat": 52.5},
        "C. St. John's, Newfoundland and Labrador": {"lat": 47.6, "lon": -52.7, "utc_offset": -3.5, "mag_lat": 56.0},
        "D. Alert, Nunavut": {"lat": 82.5, "lon": -62.3, "utc_offset": -5.0, "mag_lat": 87.0},
        "E. Thurso, Scotland": {"lat": 58.6, "lon": -3.5, "utc_offset": 0.0, "mag_lat": 61.0}
    }

    print(f"Analysis for an aurora event at {event_utc_hour:02d}:{event_utc_minute:02d} UTC with Kp={kp_index}\n")
    print("The ideal magnetic latitude for overhead aurora during a Kp=7 storm is roughly 55-60°.")
    print("--------------------------------------------------------------------------")

    # The UTC time of the event
    event_time_utc = datetime.datetime(2023, 11, 5, event_utc_hour, event_utc_minute)

    for name, data in locations.items():
        print(f"\nAnalyzing: {name}")

        # 1. Calculate Local Time
        utc_offset_hours = int(data["utc_offset"])
        utc_offset_minutes = int((data["utc_offset"] % 1) * 60)
        time_delta = datetime.timedelta(hours=utc_offset_hours, minutes=utc_offset_minutes)
        local_time = event_time_utc + time_delta
        print(f"  - Local Time: {local_time.strftime('%H:%M')}")

        # 2. Assess Darkness
        # In early November, locations A, B, C are in full darkness.
        # D (Alert) is in polar night (24h darkness).
        # E (Thurso) is in nautical twilight, with the sun rising soon, making viewing difficult.
        is_dark = True
        darkness_note = "Deep darkness, good for viewing."
        if name.startswith("E."):
            is_dark = False
            darkness_note = "In twilight before sunrise. Sky is brightening, poor for viewing."
        elif name.startswith("D."):
            darkness_note = "Polar night (24h darkness)."
        print(f"  - Viewing Conditions: {darkness_note}")

        # 3. Assess Magnetic Latitude
        mag_lat = data["mag_lat"]
        print(f"  - Magnetic Latitude: ~{mag_lat}°")
        analysis = ""
        if mag_lat > 75:
            analysis = "This is in the polar cap, north of the main auroral oval, even during a storm. Unlikely to see overhead storm aurora."
        elif 60 <= mag_lat <= 65:
            if not is_dark:
                 analysis = "Good magnetic latitude, but poor viewing time due to daylight."
            else:
                 analysis = "Good magnetic latitude, potential candidate."
        elif 55 <= mag_lat < 60:
            analysis = "Excellent magnetic latitude for an overhead display during a Kp=7 storm."
        elif 50 <= mag_lat < 55:
            analysis = "On the southern edge of the storm-expanded oval. Aurora may be visible on the horizon, but overhead is unlikely."

        print(f"  - Analysis: {analysis}")

    print("\n--------------------------------------------------------------------------")
    print("\nConclusion:")
    print(" - Portland and Madison are too far south magnetically for an overhead display.")
    print(" - Alert is too far north; the auroral oval moves south (away from Alert) during a storm.")
    print(" - Thurso, while having a good magnetic latitude, is in twilight, making aurora viewing very difficult.")
    print(" - St. John's is in deep darkness and its magnetic latitude of ~56° places it directly under the expanded auroral oval during a Kp=7 storm.")
    print("\nTherefore, St. John's is the most likely location.")

analyze_aurora_locations()