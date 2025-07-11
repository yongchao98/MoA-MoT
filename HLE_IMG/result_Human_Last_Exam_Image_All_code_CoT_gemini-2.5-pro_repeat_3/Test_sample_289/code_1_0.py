def convert_feet_to_meters():
    """
    Converts the elevation of Descanso, CA from feet to meters and prints the calculation.
    """
    # Elevation of Descanso, California in feet, according to the GNIS.
    elevation_feet = 3428

    # Conversion factor from feet to meters.
    feet_to_meters_factor = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_factor

    print(f"The elevation of Descanso, CA is {elevation_feet} feet.")
    print(f"The conversion factor from feet to meters is {feet_to_meters_factor}.")
    print("The final calculation is:")
    print(f"{elevation_feet} * {feet_to_meters_factor} = {elevation_meters}")

convert_feet_to_meters()