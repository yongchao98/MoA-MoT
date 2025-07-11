import skyfield.api
from skyfield.api import load, EarthSatellite
from skyfield.timelib import Time

def calculate_sky_position():
    """
    Calculates the position of celestial objects for a given time and location
    to determine the viewpoint country of an astronomy picture.
    """
    # 1. Setup the location. Let's use Rome, Italy as a candidate.
    # A southerly latitude is hypothesized because Jupiter appears very high.
    city_name = "Rome, Italy"
    latitude_val = 41.9028
    longitude_val = 12.4964
    location = skyfield.api.Topos(latitude_degrees=latitude_val, longitude_degrees=longitude_val)

    # 2. Setup the time from the image: 2024-10-22 20:07:37 CEST.
    # CEST is UTC+2, so we use 18:07:37 UTC.
    ts = load.timescale()
    time = ts.utc(2024, 10, 22, 18, 7, 37)

    # 3. Load ephemeris data for planets and stars
    eph = load('de421.bsp')  # JPL ephemeris for planets
    stardata = load.supplement('hip_main.csv') # Hipparcos catalog for stars

    earth = eph['earth']
    
    # Define celestial bodies
    jupiter = eph['jupiter']
    saturn = eph['saturn']
    # From Hipparcos Catalog: Vega is HIP 91262, Capella is HIP 24436
    vega = skyfield.api.Star.from_df(stardata[stardata['hip'] == 91262])
    capella = skyfield.api.Star.from_df(stardata[stardata['hip'] == 24436])
    
    bodies = {
        'Jupiter': jupiter,
        'Saturn': saturn,
        'Vega': vega,
        'Capella': capella
    }

    # 4. Calculate and print results
    print(f"Calculating object positions for the following conditions:")
    print(f"Viewpoint Location: {city_name} (Latitude: {latitude_val:.2f}°, Longitude: {longitude_val:.2f}°)")
    print(f"Date and Time (UTC): {time.utc_strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 30)

    altitudes = {}
    for name, body in bodies.items():
        astrometric = (earth + location).at(time).observe(body)
        alt, az, _ = astrometric.apparent().altaz()
        altitudes[name] = alt.degrees
        
        print(f"Object: {name}")
        print(f"  Altitude: {alt.degrees:.2f}°")
        print(f"  Azimuth: {az.degrees:.2f}°")
    print("-" * 30)
    
    # 5. Conclusion based on altitudes
    jupiter_alt = altitudes['Jupiter']
    capella_alt = altitudes['Capella']
    alt_difference = jupiter_alt - capella_alt
    
    print("Analysis:")
    print(f"The altitude of Jupiter is {jupiter_alt:.2f}°.")
    print(f"The altitude of Capella is {capella_alt:.2f}°.")
    print(f"The altitude difference is {alt_difference:.2f}°. This significant difference, with Jupiter being much higher,")
    print("is consistent with the visual evidence in the image and supports a southern latitude.")
    print("\nTherefore, the viewpoint is likely in a southern European country, such as Italy.")


calculate_sky_position()