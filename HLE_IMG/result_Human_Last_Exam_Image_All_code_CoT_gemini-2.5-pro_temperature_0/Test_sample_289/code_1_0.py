def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Descanso, California in meters.
    The elevation in feet is obtained from the Geographic Names Information System (GNIS).
    """
    # Elevation of Descanso, CA in feet (Source: GNIS)
    elevation_feet = 3451

    # Conversion factor: 1 foot = 0.3048 meters
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # Print the equation as requested
    print(f"The elevation of the unincorporated community (Descanso, CA) is {elevation_feet} feet.")
    print("To convert this to meters, we multiply by the conversion factor.")
    print(f"Calculation: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")

calculate_elevation_in_meters()