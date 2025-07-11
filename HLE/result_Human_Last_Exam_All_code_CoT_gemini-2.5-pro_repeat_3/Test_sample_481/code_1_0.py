def find_imprint_differences():
    """
    Identifies and displays the different numerical imprints on Teva brand
    fluoxetine 10mg and 20mg capsules.
    """
    # Numerical imprints for the specified capsules
    imprint_10mg = "3105"
    imprint_20mg = "3106"

    # Display the findings
    print("Based on drug identification databases:")
    print(f"The numerical imprint on the Teva 10mg fluoxetine capsule is: {imprint_10mg}")
    print(f"The numerical imprint on the Teva 20mg fluoxetine capsule is: {imprint_20mg}")
    print("\nTherefore, the different numbers imprinted on the capsules are:")
    # The final line below outputs the numbers as requested.
    print(f"{imprint_10mg} and {imprint_20mg}")

find_imprint_differences()