def find_pen_owner_family_name():
    """
    This function stores and retrieves information about the original owner
    of the lucky ballpoint pen from the Odd Taxi audio drama.
    """
    # Information from Odd Taxi audio drama 13.3 reveals the pen's origin.
    # The pen was a gift for the rakugo master, Shofukutei Donraku.
    # In Japanese naming order, the family name comes first.
    original_owner_full_name = "Shofukutei Donraku"
    
    # The family name is the first part of the full name.
    family_name = original_owner_full_name.split()[0]

    print(f"According to the Odd Taxi audio drama 13.3, the original owner of the lucky ballpoint pen was {original_owner_full_name}.")
    print(f"The family name of the original owner is: {family_name}")

find_pen_owner_family_name()