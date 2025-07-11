from skyfield.api import load, Star, Topos
from skyfield.data import hipparcos

def calculate_celestial_positions(latitude, longitude, year, month, day, hour, minute, second):
    """
    Calculates and prints the altitude and azimuth of several celestial bodies
    for a given location and time.
    """
    # Set up the timescale and the specific time
    ts = load.timescale()
    t = ts.utc(year, month, day, hour, minute, second)

    # Load ephemeris data for planets and stars
    eph = load('de421.bsp')
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)

    # Define celestial bodies
    earth = eph['earth']
    planets = {
        'Jupiter': eph['jupiter'],
        'Saturn': eph['saturn'],
        'Uranus': eph['uranus']
    }
    named_stars = {
        'Vega': Star.from_dataframe(stars.loc[91262]),
        'Capella': Star.from_dataframe(stars.loc[24608])
    }

    # Define the observer's location
    observer_location = Topos(latitude_degrees=latitude, longitude_degrees=longitude)
    observer = earth + observer_location

    # Calculate and store positions
    positions = {}
    all_bodies = {**planets, **named_stars}
    for name, body in all_bodies.items():
        astrometric = observer.at(t).observe(body)
        alt, az, _ = astrometric.apparent().altaz()
        positions[name] = {'alt': alt.degrees, 'az': az.degrees}
        
    return positions

def print_results(location_name, positions):
    """Prints the calculated positions in a formatted way."""
    print(f"--- Celestial Body Altitudes for {location_name} ---")
    highest_body = max(positions, key=lambda k: positions[k]['alt'])
    for name, pos in sorted(positions.items(), key=lambda item: item[1]['alt'], reverse=True):
        print(f"{name:>8}: Altitude = {pos['alt']:.1f}°")
    print(f"\nIn {location_name}, the highest object is {highest_body}.")
    print("-" * 45)

# --- Main execution ---
# Time from image: 2024-10-22 20:07:37 CEST (UTC+2)
utc_year, utc_month, utc_day = 2024, 10, 22
utc_hour, utc_minute, utc_second = 18, 7, 37

# Case 1: A typical European location (Berlin, Germany)
berlin_lat, berlin_lon = 52.5, 13.4
berlin_positions = calculate_celestial_positions(berlin_lat, berlin_lon, utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_second)
print_results("Berlin, Germany", berlin_positions)

# Case 2: A location in North Africa (Tripoli, Libya)
# Libya is also in a UTC+2 time zone.
tripoli_lat, tripoli_lon = 32.8, 13.2
tripoli_positions = calculate_celestial_positions(tripoli_lat, tripoli_lon, utc_year, utc_month, utc_day, utc_hour, utc_minute, utc_second)
print_results("Tripoli, Libya", tripoli_positions)
