def calculate_elevation_in_meters():
    """
    Calculates the elevation of Warner Springs, CA in meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # Elevation of Warner Springs in feet (from GNIS)
    elevation_feet = 3133

    # Conversion factor from feet to meters
    conversion_factor = 0.3048

    # Calculate elevation in meters
    elevation_meters = elevation_feet * conversion_factor
    
    # The official GNIS value in meters is 955, which is the rounded result of this calculation.
    rounded_elevation_meters = round(elevation_meters)

    print(f"The elevation of the unincorporated community (Warner Springs, CA) is calculated by converting its elevation from feet to meters.")
    print(f"Elevation in feet: {elevation_feet}")
    print(f"Conversion factor (meters per foot): {conversion_factor}")
    print(f"Calculation: {elevation_feet} * {conversion_factor} = {elevation_meters}")
    print(f"The elevation in meters, rounded to the nearest whole number, is {rounded_elevation_meters}.")

calculate_elevation_in_meters()