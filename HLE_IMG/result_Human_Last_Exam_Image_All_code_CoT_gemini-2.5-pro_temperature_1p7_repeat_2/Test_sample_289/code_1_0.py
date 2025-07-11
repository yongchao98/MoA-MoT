def convert_feet_to_meters():
    """
    Converts the elevation of Pala, CA from feet to meters and prints the calculation.
    """
    # The elevation of Pala, CA, according to the GNIS is 410 feet.
    elevation_feet = 410

    # The conversion factor from feet to meters.
    feet_to_meters_conversion_factor = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion_factor

    print(f"The official elevation of Pala, CA is {elevation_feet} feet.")
    print("To convert feet to meters, we use the conversion factor 0.3048.")
    print(f"The calculation is: {elevation_feet} * {feet_to_meters_conversion_factor} = {elevation_meters}")
    
convert_feet_to_meters()