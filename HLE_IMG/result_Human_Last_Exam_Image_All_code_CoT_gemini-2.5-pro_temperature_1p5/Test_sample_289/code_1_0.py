def calculate_elevation_in_meters():
    """
    Calculates the elevation of Descanso, California in meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # Elevation of Descanso, CA in feet, according to GNIS
    elevation_feet = 3576

    # Conversion factor from feet to meters
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # Print the calculation steps and the result
    print(f"The official elevation of the unincorporated community, Descanso, is {elevation_feet} feet.")
    print("To convert this to meters, we use the conversion factor of 1 foot = 0.3048 meters.")
    print("The final calculation is:")
    print(f"{elevation_feet} feet * {feet_to_meters_conversion} m/foot = {elevation_meters} meters")

calculate_elevation_in_meters()