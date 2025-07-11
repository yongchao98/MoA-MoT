def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Potrero, California in meters.
    The location was identified from the image. The elevation data is from the
    Geographic Names Information System (GNIS).
    """

    # The elevation of Potrero, California, in feet, according to GNIS (Feature ID: 248232).
    elevation_in_feet = 2316

    # The conversion factor from feet to meters.
    conversion_factor_feet_to_meters = 0.3048

    # Calculate the elevation in meters.
    elevation_in_meters = elevation_in_feet * conversion_factor_feet_to_meters

    print(f"The location in the image is the unincorporated community of Potrero, California.")
    print(f"According to the Geographic Names Information System, its elevation is {elevation_in_feet} feet.")
    print("To convert feet to meters, we multiply the value in feet by the conversion factor.")
    print("\nHere is the calculation:")
    print(f"{elevation_in_feet} (feet) * {conversion_factor_feet_to_meters} (m/ft) = {elevation_in_meters:.2f} (meters)")
    
    # The final answer rounded to the nearest whole number.
    final_answer = round(elevation_in_meters)
    print(f"\nThus, the elevation is approximately {final_answer} meters.")

calculate_elevation_in_meters()