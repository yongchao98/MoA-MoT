def find_elevation():
    """
    This function calculates and displays the elevation of Santa Ysabel, California in meters.
    The elevation in feet is taken from the Geographic Names Information System (GNIS).
    """
    # Elevation of Santa Ysabel, CA in feet, according to GNIS
    elevation_feet = 2995

    # Conversion factor from feet to meters
    feet_to_meters_conversion = 0.3048

    # Calculate elevation in meters
    elevation_meters = elevation_feet * feet_to_meters_conversion

    # The official GNIS value is rounded. We will use the rounded value for the final answer.
    final_elevation_meters = round(elevation_meters)

    print(f"The unincorporated community in the image is Santa Ysabel, California.")
    print(f"According to the Geographic Names Information System (GNIS), its elevation is {elevation_feet} feet.")
    print(f"To convert feet to meters, we use the conversion: 1 foot = {feet_to_meters_conversion} meters.")
    print(f"Calculation: {elevation_feet} feet * {feet_to_meters_conversion} m/ft = {elevation_meters:.3f} meters.")
    print(f"The official elevation, rounded to the nearest whole number, is {final_elevation_meters} meters.")

find_elevation()