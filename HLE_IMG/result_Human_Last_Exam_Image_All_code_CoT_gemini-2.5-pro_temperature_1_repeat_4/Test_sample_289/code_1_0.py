def get_elevation_of_descanso():
    """
    This function provides the elevation of Descanso, California,
    as found in the Geographic Names Information System (GNIS).

    The location was identified via a reverse image search, which pointed
    to Descanso, California. A subsequent search in the GNIS database
    provided the official elevation. The prompt asks for the elevation in meters.
    """
    
    # The unincorporated community identified from the image is Descanso, California.
    community_name = "Descanso, California"
    
    # According to the GNIS, the elevation of Descanso is 1042 meters.
    elevation_in_meters = 1042
    
    print(f"The unincorporated community in the image is: {community_name}")
    print(f"According to the Geographic Names Information System, the elevation as of 2014 is:")
    print(elevation_in_meters)
    print("The final answer is expressed in meters.")

get_elevation_of_descanso()