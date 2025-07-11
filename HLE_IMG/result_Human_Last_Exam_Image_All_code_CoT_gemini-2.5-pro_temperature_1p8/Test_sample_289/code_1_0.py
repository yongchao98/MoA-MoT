def find_elevation_in_meters():
    """
    This function calculates and displays the elevation of Descanso, CA in meters.
    The data is based on the Geographic Names Information System (GNIS).
    """
    # According to GNIS, the elevation of Descanso, California is 3570 feet.
    elevation_in_feet = 3570

    # The conversion factor from feet to meters.
    feet_to_meters_factor = 0.3048

    # Calculate the elevation in meters.
    elevation_in_meters = elevation_in_feet * feet_to_meters_factor
    
    # The GNIS typically provides a rounded integer value for meters.
    final_elevation = round(elevation_in_meters)

    print("The unincorporated community is Descanso, California.")
    print("According to the Geographic Names Information System (GNIS), the elevation is 3570 feet.")
    print("\nTo convert this to meters, we use the equation:")
    print(f"{elevation_in_feet} feet * {feet_to_meters_factor} meters/foot = {final_elevation} meters (rounded)")
    print(f"\nThe elevation of the community is {final_elevation} meters.")

find_elevation_in_meters()