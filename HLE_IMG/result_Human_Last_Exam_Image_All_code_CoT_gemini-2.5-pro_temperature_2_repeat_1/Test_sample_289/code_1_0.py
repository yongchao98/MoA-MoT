def calculate_elevation_in_meters():
    """
    Calculates and prints the elevation of Jamul, CA in meters based on its elevation in feet from GNIS.
    """
    # The elevation of Jamul, California in feet, according to the Geographic Names Information System (GNIS Feature ID: 272449).
    elevation_in_feet = 974

    # The conversion factor from feet to meters.
    feet_to_meters_conversion_factor = 0.3048

    # Calculate the elevation in meters.
    elevation_in_meters = elevation_in_feet * feet_to_meters_conversion_factor

    # Print the equation as requested.
    print(f"The elevation of the unincorporated community is {elevation_in_feet} feet.")
    print(f"To convert this to meters, we use the calculation:")
    print(f"{elevation_in_feet} * {feet_to_meters_conversion_factor} = {elevation_in_meters}")

calculate_elevation_in_meters()