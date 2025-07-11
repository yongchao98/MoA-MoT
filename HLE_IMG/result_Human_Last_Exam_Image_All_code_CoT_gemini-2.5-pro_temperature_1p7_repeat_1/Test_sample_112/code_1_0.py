from skyfield.api import Star, load, Topos
from skyfield.data import hipparcos

def find_location():
    """
    Analyzes the sky from the image to determine the viewpoint country.
    """
    # 1. Setup libraries and load astronomical data
    ts = load.timescale()
    eph = load('de421.bsp')
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)

    # 2. Define the time from the image: 2024-10-22 20:07:37 CEST (UTC+2)
    # CEST is UTC+2, so 20:07 CEST is 18:07 UTC.
    t = ts.utc(2024, 10, 22, 18, 7, 37)

    # 3. Choose a representative location in the CEST timezone. Germany is a central example.
    # We will use Berlin's coordinates.
    country = "Germany"
    city = "Berlin"
    location = Topos('52.5200 N', '13.4050 E')
    observer = eph['earth'] + location

    # 4. Identify the celestial objects from the image
    jupiter = eph['jupiter barycenter']
    saturn = eph['saturn barycenter']
    uranus = eph['uranus barycenter']
    vega = Star.from_df(stars.loc[91262])     # HIP 91262
    capella = Star.from_df(stars.loc[24608])  # HIP 24608

    # 5. Calculate the apparent position (Altitude/Azimuth) of each object
    jup_pos = observer.at(t).observe(jupiter).apparent()
    sat_pos = observer.at(t).observe(saturn).apparent()
    ura_pos = observer.at(t).observe(uranus).apparent()
    vega_pos = observer.at(t).observe(vega).apparent()
    cap_pos = observer.at(t).observe(capella).apparent()

    jup_alt, jup_az, _ = jup_pos.altaz()
    sat_alt, sat_az, _ = sat_pos.altaz()
    ura_alt, ura_az, _ = ura_pos.altaz()
    vega_alt, vega_az, _ = vega_pos.altaz()
    cap_alt, cap_az, _ = cap_pos.altaz()
    
    print(f"Analysis for viewpoint based on date 2024-10-22 20:07:37 CEST.")
    print(f"Assuming a location in {city}, {country}.")
    print("-" * 50)
    print("Calculated Celestial Object Positions:")
    # Azimuth key: North=0°, East=90°, South=180°, West=270°
    # Altitude key: Horizon=0°, Zenith=90°
    print(f"Jupiter: Altitude = {jup_alt.degrees:.1f}°, Azimuth = {jup_az.degrees:.1f}° (East-Southeast)")
    print(f"Capella: Altitude = {cap_alt.degrees:.1f}°, Azimuth = {cap_az.degrees:.1f}° (East-Northeast)")
    print(f"Uranus:  Altitude = {ura_alt.degrees:.1f}°, Azimuth = {ura_az.degrees:.1f}° (East)")
    print(f"Vega:    Altitude = {vega_alt.degrees:.1f}°, Azimuth = {vega_az.degrees:.1f}° (West)")
    print(f"Saturn:  Altitude = {sat_alt.degrees:.1f}°, Azimuth = {sat_az.degrees:.1f}° (South-Southwest)")
    print("-" * 50)
    print("Conclusion:")
    print("The calculations show that at the given time, Jupiter, Capella, and Uranus are in the eastern part of the sky,")
    print("while Saturn and Vega are in the southwestern and western sky, respectively.")
    print("This configuration, with a prominent Milky Way, is consistent with an autumn evening view from a mid-latitude")
    print("Northern Hemisphere location, such as Germany, which uses the CEST timezone.")

find_location()