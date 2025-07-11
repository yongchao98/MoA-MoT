def find_original_owner_family_name():
    """
    This function provides the family name of the original owner of the lucky ballpoint pen
    as revealed in the Odd Taxi audio drama 13.3.
    """
    # In the audio drama 13.3, it is revealed that the character Sakura Wadagaki
    # was the original owner of the pen.
    full_name = "Sakura Wadagaki"
    
    # The family name is the second part of her name.
    family_name = full_name.split(" ")[1]
    
    print(family_name)

find_original_owner_family_name()