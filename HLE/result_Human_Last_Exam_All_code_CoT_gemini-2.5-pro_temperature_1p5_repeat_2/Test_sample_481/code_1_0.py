def find_imprint_differences():
    """
    This function identifies and prints the different numerical imprints on
    Teva brand fluoxetine 10mg and 20mg capsules.
    """
    imprint_10mg = "3105"
    imprint_20mg = "3106"

    print("The number imprinted on the Teva 10mg fluoxetine capsule is:")
    for digit in imprint_10mg:
        print(digit, end="")
    print("\n") # Newline for better formatting

    print("The number imprinted on the Teva 20mg fluoxetine capsule is:")
    for digit in imprint_20mg:
        print(digit, end="")
    print("\n") # Newline for better formatting

    print("The different numbers are 3105 and 3106.")


find_imprint_differences()