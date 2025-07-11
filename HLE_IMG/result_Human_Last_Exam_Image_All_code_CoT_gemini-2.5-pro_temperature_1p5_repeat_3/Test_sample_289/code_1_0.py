def convert_feet_to_meters():
    """
    Converts the elevation of Jamul, CA from feet to meters.
    The official elevation is 965 feet according to GNIS.
    The conversion factor is 1 foot = 0.3048 meters.
    """
    elevation_feet = 965
    feet_to_meters_conversion_factor = 0.3048
    
    # Calculate the elevation in meters
    elevation_meters = elevation_feet * feet_to_meters_conversion_factor
    
    # Print the equation and the result
    print(f"The elevation of Jamul is {elevation_feet} feet.")
    print(f"To convert to meters, we use the formula: Elevation_in_meters = Elevation_in_feet * {feet_to_meters_conversion_factor}")
    print(f"Calculation: {elevation_feet} * {feet_to_meters_conversion_factor} = {elevation_meters:.3f}")
    print(f"The elevation of Jamul is approximately {round(elevation_meters)} meters.")

convert_feet_to_meters()