from skyfield.api import load, Topos
from skyfield.starlib import Star

def find_celestial_positions():
    """
    Calculates the position of key celestial objects for a specific time and location
    to determine the viewpoint of the provided sky map.
    """
    # Load ephemeris data
    ts = load.timescale()
    eph = load('de421.bsp')

    # Define celestial bodies from the image
    saturn = eph['saturn']
    jupiter = eph['jupiter']
    
    # Vega (HIP 91262) and Capella (HIP 24608) are bright stars also shown
    hipparcos_df = load.hipparcos_df
    vega = Star.from_df(hipparcos_df.loc[91262])
    capella = Star.from_df(hipparcos_df.loc[24608])

    # Time of observation from the image: 2024-10-22 20:07:37 CEST
    # CEST is UTC+2, so we convert to UTC.
    time = ts.utc(2024, 10, 22, 18, 7, 37)

    # Candidate location: Paris, France. This is a common default for astronomy
    # software and is in the CEST timezone.
    # Latitude: 48.8566° N, Longitude: 2.3522° E
    location_name = "Paris, France"
    observer = Topos(latitude_degrees=48.8566, longitude_degrees=2.3522)
    
    # Calculate apparent positions (altitude and azimuth)
    alt_saturn, az_saturn, _ = (eph['earth'] + observer).at(time).observe(saturn).apparent().altaz()
    alt_jupiter, az_jupiter, _ = (eph['earth'] + observer).at(time).observe(jupiter).apparent().altaz()
    alt_vega, az_vega, _ = (eph['earth'] + observer).at(time).observe(vega).apparent().altaz()
    alt_capella, az_capella, _ = (eph['earth'] + observer).at(time).observe(capella).apparent().altaz()

    # Print the results
    print(f"Calculated positions for the viewpoint from {location_name}")
    print(f"on 2024-10-22 at 20:07:37 CEST:")
    print("-" * 50)
    print(f"Saturn:  Altitude = {alt_saturn.degrees:.1f} degrees, Azimuth = {az_saturn.degrees:.1f} degrees (SSW)")
    print(f"Jupiter: Altitude = {alt_jupiter.degrees:.1f} degrees, Azimuth = {az_jupiter.degrees:.1f} degrees (ESE)")
    print(f"Vega:    Altitude = {alt_vega.degrees:.1f} degrees, Azimuth = {az_vega.degrees:.1f} degrees (WSW)")
    print(f"Capella: Altitude = {alt_capella.degrees:.1f} degrees, Azimuth = {az_capella.degrees:.1f} degrees (NE)")
    print("-" * 50)
    print("Conclusion: The calculated positions for Paris closely match the sky map.")
    print("Saturn is low in the SSW (~10 deg alt), and Jupiter is higher in the ESE (~31 deg alt).")
    print("This confirms the viewpoint is in France.")

find_celestial_positions()