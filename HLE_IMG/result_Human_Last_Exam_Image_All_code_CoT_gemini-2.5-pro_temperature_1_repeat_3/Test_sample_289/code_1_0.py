def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Descanso, CA in meters.
    The official elevation in feet is obtained from the Geographic Names Information System (GNIS).
    """
    # Elevation of Descanso, CA in feet according to GNIS.
    elevation_feet = 3442

    # Conversion factor from feet to meters.
    feet_to_meters_conversion = 0.3048

    # Calculate elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # The GNIS also directly provides the elevation in meters, which is 1049.
    # We will use the directly provided value for higher accuracy as it's the official record.
    official_elevation_meters = 1049

    print(f"The official elevation of Descanso, California is {official_elevation_meters} meters.")

calculate_elevation_in_meters()