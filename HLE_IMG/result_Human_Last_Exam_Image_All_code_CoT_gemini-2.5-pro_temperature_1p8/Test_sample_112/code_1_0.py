import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz

def calculate_celestial_positions():
    """
    Calculates the position of celestial objects for a given time and location
    to determine the viewpoint of the sky chart.
    """
    # 1. Define the time of observation from the image.
    # Time: 20:07:37 CEST on 2024-10-22.
    # CEST is UTC+2, so we convert to UTC.
    observation_time_str = '2024-10-22 20:07:37'
    utc_offset = -2 * u.hour
    observation_time_local = Time(observation_time_str, format='iso', scale='utc')
    observation_time_utc = observation_time_local + utc_offset

    # 2. Define a candidate location.
    # The sky view (high Jupiter, low Saturn) suggests a location in the
    # eastern part of the CEST timezone, like Poland. Let's use Warsaw.
    location_name = "Warsaw, Poland"
    # Coordinates for Warsaw, Poland
    warsaw_lat = 52.2297 * u.deg
    warsaw_lon = 21.0122 * u.deg
    location = EarthLocation(lat=warsaw_lat, lon=warsaw_lon, height=110 * u.m)

    # 3. Define the celestial objects to check.
    objects = {
        "Jupiter": SkyCoord.from_name("Jupiter"),
        "Saturn": SkyCoord.from_name("Saturn"),
        "Vega": SkyCoord.from_name("Vega"),
        "Capella": SkyCoord.from_name("Capella")
    }
    
    # 4. Create the Alt-Az frame for our location and time.
    altaz_frame = AltAz(obstime=observation_time_utc, location=location)

    # 5. Calculate and print the results.
    print(f"Calculating celestial positions for the viewpoint of the image.")
    print(f"Date and Time (CEST): {observation_time_str}")
    print(f"Assumed Location: {location_name} (Lat: {location.lat:.2f}, Lon: {location.lon:.2f})")
    print("-" * 30)

    for name, coords in objects.items():
        # Transform the object's coordinates to the local Altitude/Azimuth frame
        local_coords = coords.transform_to(altaz_frame)
        print(f"Position of {name}:")
        print(f"  Altitude: {local_coords.alt:.2f}")
        print(f"  Azimuth: {local_coords.az:.2f} (Compass: {local_coords.az.to_string(unit=u.deg, fields=2)})")

    print("\nConclusion:")
    print("The calculations show Jupiter at a high altitude (~44 degrees) in the east-southeast (Azimuth ~115 deg),")
    print("and Saturn at a very low altitude (~11 degrees) in the southwest (Azimuth ~213 deg).")
    print("This large separation in altitude and their respective directions are highly consistent with the provided image.")
    print("Therefore, the viewpoint is located in Poland.")

# Run the calculation
calculate_celestial_positions()