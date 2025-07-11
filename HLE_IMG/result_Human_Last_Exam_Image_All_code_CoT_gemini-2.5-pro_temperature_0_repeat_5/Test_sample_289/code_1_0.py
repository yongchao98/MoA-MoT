def convert_elevation_to_meters():
    """
    This function calculates the elevation of Descanso, CA in meters.
    The location was identified from the image. The elevation data is from the
    Geographic Names Information System (GNIS).
    """
    # According to the GNIS, the elevation of Descanso, CA is 3514 feet.
    elevation_feet = 3514

    # The conversion factor from feet to meters is 0.3048.
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # Print the steps of the conversion.
    print(f"The elevation of the unincorporated community (Descanso, CA) is {elevation_feet} feet.")
    print("The conversion to meters is calculated by multiplying the elevation in feet by the conversion factor.")
    print(f"Calculation: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")
    print(f"The elevation in meters is {round(elevation_meters)}.")

convert_elevation_to_meters()