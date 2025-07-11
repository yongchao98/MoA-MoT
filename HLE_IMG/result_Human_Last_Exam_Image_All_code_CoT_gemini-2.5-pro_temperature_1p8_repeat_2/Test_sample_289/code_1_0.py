def calculate_elevation_in_meters():
    """
    This function converts the elevation of Santa Ysabel, California, from feet to meters.
    The elevation in feet is sourced from the Geographic Names Information System (GNIS).
    """
    # The official elevation of Santa Ysabel, CA in feet according to GNIS.
    elevation_feet = 2986

    # The conversion factor from feet to meters.
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion

    print(f"The elevation of the unincorporated community, Santa Ysabel, is {elevation_feet} feet.")
    print(f"To convert this to meters, we multiply by the conversion factor {feet_to_meters_conversion}.")
    print(f"The calculation is: {elevation_feet} * {feet_to_meters_conversion} = {elevation_meters}")
    print(f"The elevation in meters is {elevation_meters:.2f}.")

calculate_elevation_in_meters()