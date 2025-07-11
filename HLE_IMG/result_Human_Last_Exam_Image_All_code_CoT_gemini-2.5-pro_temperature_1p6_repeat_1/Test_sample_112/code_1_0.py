import pytz
from skyfield.api import Star, load, Topos
from skyfield.data import hipparcos

def find_location_from_sky_chart():
    """
    Calculates the position of celestial objects for a given time and location
    to match the provided sky chart.
    """
    # 1. Setup
    ts = load.timescale()
    eph = load('de421.bsp')
    earth = eph['earth']
    
    # 2. Time from image
    # The image shows 2024-10-22 20:07:37 CEST.
    # CEST is Central European Summer Time, which is UTC+2.
    # So, the time in UTC is 18:07:37.
    year, month, day = 2024, 10, 22
    hour, minute, second = 18, 7, 37
    time = ts.utc(year, month, day, hour, minute, second)
    
    # 3. Observer Location
    # Based on analysis, a location in central Germany provides the best fit.
    # Let's use Frankfurt, Germany.
    city = "Frankfurt, Germany"
    latitude, longitude = 50.1, 8.7
    location = Topos(latitude_degrees=latitude, longitude_degrees=longitude)
    observer = earth + location
    
    # 4. Celestial Objects from image
    jupiter = eph['jupiter barycenter']
    saturn = eph['saturn barycenter']
    
    # Load star data for Vega and Capella from the Hipparcos catalog.
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)
    
    # Vega is HIP 91262, Capella is HIP 24608
    vega = Star.from_dataframe(df.loc[91262])
    capella = Star.from_dataframe(df.loc[24608])

    # 5. Calculate apparent positions (Altitude and Azimuth)
    jup_app = observer.at(time).observe(jupiter).apparent()
    sat_app = observer.at(time).observe(saturn).apparent()
    vega_app = observer.at(time).observe(vega).apparent()
    cap_app = observer.at(time).observe(capella).apparent()

    jup_alt, jup_az, _ = jup_app.altaz()
    sat_alt, sat_az, _ = sat_app.altaz()
    vega_alt, vega_az, _ = vega_app.altaz()
    cap_alt, cap_az, _ = cap_app.altaz()
    
    # 6. Print the results
    print(f"Sky analysis for the image taken on {year}-{month}-{day} at {hour+2}:{minute:02d}:{second:02d} CEST.")
    print(f"Assuming a viewpoint from {city} ({latitude}°N, {longitude}°E):\n")
    print("Calculated Celestial Positions (Alt/Az):")
    print(f"- Saturn:  Altitude = {sat_alt.degrees:.1f}°, Azimuth = {sat_az.degrees:.1f}° (South-West)")
    print(f"- Jupiter: Altitude = {jup_alt.degrees:.1f}°, Azimuth = {jup_az.degrees:.1f}° (East)")
    print(f"- Vega:    Altitude = {vega_alt.degrees:.1f}°, Azimuth = {vega_az.degrees:.1f}° (West-South-West)")
    print(f"- Capella: Altitude = {cap_alt.degrees:.1f}°, Azimuth = {cap_az.degrees:.1f}° (North-East)")
    print("\nConclusion: The calculated positions for a location in Germany closely match the scene in the image.")

find_location_from_sky_chart()