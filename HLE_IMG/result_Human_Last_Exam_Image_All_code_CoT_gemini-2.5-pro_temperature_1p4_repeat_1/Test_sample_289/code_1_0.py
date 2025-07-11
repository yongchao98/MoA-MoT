def calculate_elevation_in_meters():
    """
    This function calculates the elevation of Descanso, CA in meters.
    The elevation in feet is sourced from the USGS Geographic Names Information System.
    """
    # According to the USGS GNIS, the elevation of Descanso, CA is 3524 feet.
    elevation_feet = 3524

    # The conversion factor from feet to meters is 0.3048.
    feet_to_meters_conversion = 0.3048

    # Calculate the elevation in meters.
    elevation_meters = elevation_feet * feet_to_meters_conversion

    print("The unincorporated community in the image is Descanso, California.")
    print(f"According to the Geographic Names Information System, its elevation is {elevation_feet} feet.")
    print(f"To convert feet to meters, we multiply by the conversion factor {feet_to_meters_conversion}.")
    print("\nHere is the calculation:")
    print(f"{elevation_feet} feet * {feet_to_meters_conversion} meters/foot = {elevation_meters} meters")
    
    # Returning the final value for the answer block, but the script prints the full context.
    return elevation_meters

if __name__ == "__main__":
    calculate_elevation_in_meters()