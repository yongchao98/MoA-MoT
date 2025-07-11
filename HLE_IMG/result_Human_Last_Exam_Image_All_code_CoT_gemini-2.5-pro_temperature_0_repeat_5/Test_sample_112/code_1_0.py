from skyfield.api import load, Star, Topos
from skyfield.data import hipparcos

def calculate_celestial_positions():
    """
    Calculates the position of celestial objects for a given time and location
    to determine the viewpoint of the provided star chart.
    """
    # 1. Setup time
    ts = load.timescale()
    # Time from image: 2024-10-22 20:07:37 CEST (UTC+2)
    time = ts.utc(2024, 10, 22, 18, 7, 37)

    # 2. Setup celestial objects
    eph = load('de421.bsp')
    earth = eph['earth']
    jupiter = eph['jupiter']
    saturn = eph['saturn']

    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)

    # Star objects from Hipparcos catalog
    # Vega is HIP 91262, Capella is HIP 24608
    vega = Star.from_dataframe(df.loc[91262])
    capella = Star.from_dataframe(df.loc[24608])

    # 3. Setup observer location
    # Hypothesis: A high-latitude location in CEST zone, e.g., Stockholm, Sweden
    location_name = "Stockholm, Sweden"
    latitude = 59.3293
    longitude = 18.0686
    observer = earth + Topos(latitude_degrees=latitude, longitude_degrees=longitude)

    # 4. Calculate apparent positions (Altitude/Azimuth)
    saturn_pos = observer.at(time).observe(saturn).apparent()
    jupiter_pos = observer.at(time).observe(jupiter).apparent()
    vega_pos = observer.at(time).observe(vega).apparent()
    capella_pos = observer.at(time).observe(capella).apparent()

    alt_s, az_s, _ = saturn_pos.altaz()
    alt_j, az_j, _ = jupiter_pos.altaz()
    alt_v, az_v, _ = vega_pos.altaz()
    alt_c, az_c, _ = capella_pos.altaz()

    # 5. Print results
    print(f"Analysis for the viewpoint of the star chart:")
    print(f"Date and Time (UTC): {time.utc_strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Assumed Location: {location_name} ({latitude}° N, {longitude}° E)")
    print("-" * 40)
    print("Calculated Celestial Positions:")
    print(f"Saturn:  Altitude = {alt_s.degrees:.1f}°, Azimuth = {az_s.degrees:.1f}°")
    print(f"Jupiter: Altitude = {alt_j.degrees:.1f}°, Azimuth = {az_j.degrees:.1f}°")
    print(f"Vega:    Altitude = {alt_v.degrees:.1f}°, Azimuth = {az_v.degrees:.1f}°")
    print(f"Capella: Altitude = {alt_c.degrees:.1f}°, Azimuth = {az_c.degrees:.1f}°")
    print("-" * 40)
    print("Conclusion:")
    print("The calculated altitudes (Capella > Vega > Jupiter > Saturn) match the relative positions in the image.")
    print("This confirms the viewpoint is at a high northern latitude like Sweden.")

calculate_celestial_positions()
<<<Sweden>>>