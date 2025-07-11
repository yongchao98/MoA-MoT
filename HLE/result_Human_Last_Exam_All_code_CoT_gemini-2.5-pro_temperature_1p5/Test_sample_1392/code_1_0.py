import math

def get_local_time(utc_hour, utc_minute, offset):
    """Calculates local time given a UTC time and a timezone offset."""
    # Total minutes from midnight UTC
    total_utc_minutes = utc_hour * 60 + utc_minute
    # Total minutes from midnight local time
    local_minutes = total_utc_minutes + offset * 60
    # Handle day wrap-around
    local_minutes %= 1440
    
    local_hour = local_minutes // 60
    local_minute = local_minutes % 60
    
    return f"{int(local_hour):02d}:{int(local_minute):02d}"

def analyze_locations():
    """
    Analyzes locations to find the most likely to see an overhead aurora
    during a Kp=7 event at 06:30 UTC.
    """
    # Location data: Name, Timezone offset from UTC, approx. Geomagnetic Latitude
    locations = {
        'A': {'name': "Portland, Oregon", 'tz_offset': -8, 'mag_lat': 51.5},
        'B': {'name': "Madison, Wisconsin", 'tz_offset': -6, 'mag_lat': 53.5},
        'C': {'name': "St. John's, Newfoundland", 'tz_offset': -3.5, 'mag_lat': 58.0},
        'D': {'name': "Alert, Nunavut", 'tz_offset': -5, 'mag_lat': 86.0},
        'E': {'name': "Thurso, Scotland", 'tz_offset': 0, 'mag_lat': 58.5},
    }
    
    utc_hour = 6
    utc_minute = 30

    print(f"Analysis for an overhead aurora at {utc_hour:02d}:{utc_minute:02d} UTC during a Kp=7 event.\n")
    print("An overhead aurora during a Kp=7 storm is most likely at a geomagnetic latitude of ~50-55° and near local midnight.\n")

    for key, loc in locations.items():
        local_time_str = get_local_time(utc_hour, utc_minute, loc['tz_offset'])
        mag_lat = loc['mag_lat']

        print(f"--- Analyzing Option {key}: {loc['name']} ---")
        print(f"Local Time: {local_time_str} | Geomagnetic Latitude: ~{mag_lat}°")

        # Analysis logic
        if mag_lat > 75:
            print("Verdict: Unlikely. This location is deep within the polar cap, north of the main auroral activity during a storm.")
        elif 52 <= mag_lat <= 55:
            if 22 <= float(local_time_str.split(':')[0]) or float(local_time_str.split(':')[0]) <= 2:
                print("Verdict: Most Likely. The location is at an ideal geomagnetic latitude and is near local midnight, the time of peak auroral intensity.")
            else:
                print("Verdict: Possible. The latitude is ideal, but the timing is not optimal for peak brightness.")
        elif mag_lat > 55:
            print("Verdict: Less Likely. While this latitude sees aurora, it may be slightly north of the brightest part of the storm. The main activity would be to the south.")
        else: # mag_lat < 52
            print("Verdict: Unlikely. The location is likely too far south. The aurora might be visible on the horizon, but not directly overhead.")
        print("-" * (31 + len(loc['name'])))
        print()

analyze_locations()