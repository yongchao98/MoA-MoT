def find_imprint_differences():
    """
    This function identifies and prints the different numerical imprints
    on Teva brand fluoxetine 10mg and 20mg capsules.
    """
    
    # Numerical imprint on the 10mg capsule
    imprint_10mg = 3105
    
    # Numerical imprint on the 20mg capsule
    imprint_20mg = 3106
    
    print(f"The number imprinted on the Teva brand fluoxetine 10mg capsule is: {imprint_10mg}")
    print(f"The number imprinted on the Teva brand fluoxetine 20mg capsule is: {imprint_20mg}")
    
    print("\nThe final equation showing the different numbers is:")
    print(f"{imprint_10mg} != {imprint_20mg}")

find_imprint_differences()