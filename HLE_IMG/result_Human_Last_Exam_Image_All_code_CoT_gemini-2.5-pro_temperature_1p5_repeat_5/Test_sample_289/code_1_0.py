def convert_feet_to_meters():
    """
    Converts the elevation of Pala, California from feet to meters
    and prints the calculation.
    """
    # Elevation of Pala, CA in feet, from the Geographic Names Information System (GNIS).
    elevation_in_feet = 410

    # The standard conversion factor for feet to meters.
    conversion_factor = 0.3048

    # Calculate the elevation in meters.
    elevation_in_meters = elevation_in_feet * conversion_factor

    # Print the equation with each number.
    print(f"The elevation of the unincorporated community of Pala, CA is {elevation_in_feet} feet.")
    print("The final calculation to convert feet to meters is:")
    print(f"{elevation_in_feet} * {conversion_factor} = {elevation_in_meters}")

convert_feet_to_meters()