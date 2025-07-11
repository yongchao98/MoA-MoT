def calculate_elevation_in_meters():
    """
    Calculates the elevation of Pine Valley, CA in meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # According to GNIS (ID: 247854), the elevation of Pine Valley, CA is 3737 feet.
    elevation_in_feet = 3737
    
    # The conversion factor from feet to meters is 0.3048.
    feet_to_meters_conversion_factor = 0.3048
    
    # Perform the conversion.
    elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor
    
    # Print the equation showing all the numbers.
    print(f"The elevation of the unincorporated community (Pine Valley, CA) is {elevation_in_feet} feet.")
    print("To convert from feet to meters, we use the conversion factor of 0.3048.")
    print(f"Calculation: {elevation_in_feet} * {feet_to_meters_conversion_factor} = {elevation_in_meters}")
    
    # The final answer should be an integer, so we round the result.
    final_elevation = round(elevation_in_meters)
    print(f"\nThe elevation in meters, rounded to the nearest whole number, is {final_elevation}.")

calculate_elevation_in_meters()