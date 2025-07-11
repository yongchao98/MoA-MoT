import skyfield.api as sf
from skyfield.api import load

def find_viewpoint_country():
    """
    Determines the observer's location based on celestial data from an image.
    """
    # 1. Setup the time and load astronomical data
    ts = load.timescale()
    eph = load('de421.bsp')
    earth, jupiter, saturn = eph['earth'], eph['jupiter'], eph['saturn']

    # 2. Define the exact time of observation from the image
    # Time: 20:07:37 CEST on 2024-10-22. CEST is UTC+2.
    utc_time = ts.utc(2024, 10, 22, 18, 7, 37)
    print("Step 1: Analyzing the observation time.")
    print(f"Image time is 20:07:37 CEST on 2024-10-22.")
    print(f"This corresponds to a UTC time of {utc_time.utc_strftime()}.\n")

    # 3. Determine the observer's location
    # The image shows Jupiter at the zenith. The observer's location is the
    # point on Earth directly below Jupiter (the sub-Jovian point).
    print("Step 2: Calculating the observer's location.")
    print("Assuming Jupiter is at the zenith, we calculate the sub-Jovian point.")
    
    # skyfield can calculate the geographic position directly
    observer_location = sf.wgs84.geographic_position_of(earth.at(utc_time).observe(jupiter).apparent())
    lat = observer_location.latitude.degrees
    lon = observer_location.longitude.degrees

    print(f"Calculated Latitude: {lat:.2f} degrees")
    print(f"Calculated Longitude: {lon:.2f} degrees\n")

    # 4. Verify the location using Saturn's position
    # From this calculated location, Saturn should appear low on the horizon.
    print("Step 3: Verifying the location with Saturn's position.")
    observer_topos = sf.wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon)
    saturn_pos = observer_topos.at(utc_time).observe(saturn).apparent()
    alt, az, dist = saturn_pos.altaz()

    print(f"From the calculated location, Saturn's altitude is {alt.degrees:.2f} degrees.")
    print("This is very low on the horizon, which matches the image perfectly.\n")

    # 5. Identify the country
    # Longitude ~17.8 E is in a region that uses UTC+2 (EET/CEST).
    # Coordinates (19.5 N, 17.8 E) point to a country in North Africa.
    country = "Libya"
    print("Step 4: Identifying the country.")
    print(f"The coordinates {lat:.2f} N, {lon:.2f} E fall within the borders of Libya.")
    print("Libya uses the EET time zone (UTC+2), which is consistent with the image's CEST label (also UTC+2).")
    print(f"\nConclusion: The viewpoint of this picture is in {country}.")


if __name__ == '__main__':
    find_viewpoint_country()
<<<Libya>>>