def dms_to_dd(degrees, minutes, seconds, direction):
    """
    Converts coordinates from DMS (Degrees, Minutes, Seconds) format
    to DD (Decimal Degrees) format.
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60)
    if direction in ['S', 'W']:
        dd *= -1
    return dd

def get_site_information():
    """
    Provides information about the site at the given coordinates.
    """
    # Coordinates from the user prompt
    lat_degrees, lat_minutes, lat_seconds, lat_direction = 29, 6, 18.75, 'N'
    lon_degrees, lon_minutes, lon_seconds, lon_direction = 103, 47, 50.28, 'W'

    # Convert for potential use, and for display
    lat_dd = dms_to_dd(lat_degrees, lat_minutes, lat_seconds, lat_direction)
    lon_dd = dms_to_dd(lon_degrees, lon_minutes, lon_seconds, lon_direction)

    print(f"Analyzing location for coordinates:")
    print(f"Latitude: {lat_degrees}° {lat_minutes}' {lat_seconds}'' {lat_direction}")
    print(f"Longitude: {lon_degrees}° {lon_minutes}' {lon_seconds}'' {lon_direction}")
    print(f"(Decimal: {lat_dd}, {lon_dd})")
    print("-" * 40)

    # Information gathered about the location
    site_name = "The alleged 'Coyame UFO Crash Site'"
    location_info = "The location is in the Chihuahua Desert, near Coyame, Mexico."
    structure_info = (
        "The square feature in the image is widely believed by ufologists to be a concrete helipad or "
        "the foundation of a structure built during a covert military operation. According to a "
        "prominent conspiracy theory from 1974, a UFO and a small plane collided in this area. "
        "The story claims a secret US team was sent to recover the wreckage."
    )
    landmark_status = (
        "It is NOT an official historic landmark recognized by any government or historical preservation "
        "agency. Its notoriety comes from modern folklore and UFO conspiracy theories. The structure itself "
        "could have a more mundane origin, such as being part of an abandoned ranch or a survey site."
    )

    print(f"What is this structure? \n{structure_info}\n")
    print(f"Is it a historic landmark? \n{landmark_status}")

if __name__ == '__main__':
    get_site_information()