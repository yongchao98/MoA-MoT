def calculate_elevation_in_meters():
    """
    Converts the elevation of Descanso, CA from feet to meters.
    According to the Geographic Names Information System (GNIS),
    the elevation is 3524 feet.
    """
    elevation_feet = 3524
    conversion_factor = 0.3048  # 1 foot = 0.3048 meters

    elevation_meters = elevation_feet * conversion_factor

    print(f"The elevation of Descanso, CA is {elevation_feet} feet.")
    print("To convert feet to meters, we use the conversion factor: 1 foot = 0.3048 meters.")
    print("\nCalculation:")
    print(f"{elevation_feet} feet * {conversion_factor} m/foot = {elevation_meters:.4f} meters")

calculate_elevation_in_meters()