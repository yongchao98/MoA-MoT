def calculate_elevation_in_meters():
    """
    This script calculates the elevation of Aguanga, CA in meters.
    The elevation in feet is obtained from the Geographic Names Information System (GNIS).
    """
    # The elevation of Aguanga, California, in feet, according to GNIS (Feature ID: 238350)
    elevation_in_feet = 1955

    # The conversion factor from feet to meters
    feet_to_meters_conversion_factor = 0.3048

    # Calculate the elevation in meters
    elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor

    # Print the equation and the result
    print(f"The elevation of the unincorporated community (Aguanga, CA) is {elevation_in_feet} feet.")
    print("The conversion from feet to meters is calculated as follows:")
    print(f"{elevation_in_feet} feet * {feet_to_meters_conversion_factor} meters/foot = {elevation_in_meters} meters")

calculate_elevation_in_meters()