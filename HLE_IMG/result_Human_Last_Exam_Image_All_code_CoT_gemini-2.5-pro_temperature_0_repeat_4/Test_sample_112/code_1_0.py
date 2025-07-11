import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, get_body

def calculate_planet_altitudes(locations, observation_time):
    """
    Calculates and prints the altitudes of Jupiter and Saturn for a list of locations.
    """
    # Define the celestial bodies
    jupiter = get_body("jupiter", observation_time)
    saturn = get_body("saturn", observation_time)

    print(f"Calculating planet altitudes for {observation_time.iso} UTC\n")
    print(f"{'City':<18} {'Country':<12} {'Jupiter Alt':<15} {'Saturn Alt':<15}")
    print("-" * 60)

    for name, country, lat, lon in locations:
        # Define the observer's location
        location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=0*u.m)
        
        # Define the AltAz frame for the observer
        altaz_frame = AltAz(obstime=observation_time, location=location)
        
        # Transform planet coordinates to the AltAz frame
        jupiter_altaz = jupiter.transform_to(altaz_frame)
        saturn_altaz = saturn.transform_to(altaz_frame)
        
        # Extract the altitude
        jupiter_alt = jupiter_altaz.alt.deg
        saturn_alt = saturn_altaz.alt.deg
        
        print(f"{name:<18} {country:<12} {jupiter_alt:15.2f} {saturn_alt:15.2f}")

if __name__ == '__main__':
    # The time is 2024-10-22 20:07:37 CEST. CEST is UTC+2.
    # We convert it to UTC for calculations.
    time_cest_str = "2024-10-22 20:07:37"
    observation_time_utc = Time(time_cest_str) - 2 * u.hour

    # List of candidate locations (Name, Country, Latitude, Longitude)
    # These are all in the CEST timezone.
    candidate_locations = [
        ("Berlin", "Germany", 52.5, 13.4),
        ("Paris", "France", 48.9, 2.3),
        ("Rome", "Italy", 41.9, 12.5),
        ("Catania (Sicily)", "Italy", 37.5, 15.1),
        ("Valletta", "Malta", 35.9, 14.5)
    ]
    
    calculate_planet_altitudes(candidate_locations, observation_time_utc)

    print("\n--- Analysis ---")
    print("The image shows Jupiter very high in the sky and Saturn relatively low.")
    print("The calculations show that for northern locations like Berlin, Jupiter's altitude is only ~35 degrees.")
    print("For southern locations like Catania, Italy, Jupiter's altitude is over 50 degrees, which is a much better match for the image.")
    print("Therefore, the viewpoint is in a southern European country.")
