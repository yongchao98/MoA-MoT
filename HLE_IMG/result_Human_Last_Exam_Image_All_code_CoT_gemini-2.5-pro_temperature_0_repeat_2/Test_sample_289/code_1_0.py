def calculate_elevation_in_meters():
    """
    This script calculates the elevation of Descanso, CA in meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # Elevation of Descanso, CA in feet, according to GNIS (Feature ID 241131)
    elevation_feet = 3520

    # Conversion factor from feet to meters
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters
    elevation_meters = elevation_feet * feet_to_meters_conversion

    print(f"The elevation of the unincorporated community (Descanso, CA) is {elevation_feet} feet.")
    print(f"To convert feet to meters, we use the calculation: {elevation_feet} * {feet_to_meters_conversion}")
    print(f"The elevation in meters is: {elevation_meters:.2f}")
    # The final answer should be an integer.
    print(f"Rounding to the nearest whole number, the elevation is {round(elevation_meters)} meters.")

calculate_elevation_in_meters()