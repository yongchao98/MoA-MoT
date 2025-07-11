def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Pala, California in meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # According to the GNIS, the elevation of Pala, California is 404 feet.
    elevation_feet = 404

    # The conversion factor from feet to meters.
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # Print the full equation and the final result.
    print(f"The elevation of the unincorporated community (Pala, CA) is {elevation_feet} feet.")
    print("To convert this to meters, we multiply by the conversion factor 0.3048.")
    print(f"The calculation is: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")

calculate_elevation_in_meters()