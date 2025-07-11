def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Descanso, California in meters.
    The elevation in feet is sourced from the USGS Geographic Names Information System (GNIS).
    """
    # According to the GNIS, the elevation of Descanso, CA (Feature ID: 241113) is 3,451 feet.
    elevation_in_feet = 3451

    # The conversion factor from feet to meters.
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters.
    elevation_in_meters = elevation_in_feet * feet_to_meters_conversion

    print(f"The elevation of the unincorporated community, Descanso, is {elevation_in_feet} feet.")
    print(f"To convert this to meters, we use the conversion factor of {feet_to_meters_conversion} meters per foot.")
    print(f"The calculation is: {elevation_in_feet} * {feet_to_meters_conversion} = {elevation_in_meters}")
    print(f"Therefore, the elevation is {elevation_in_meters:.2f} meters.")

calculate_elevation_in_meters()