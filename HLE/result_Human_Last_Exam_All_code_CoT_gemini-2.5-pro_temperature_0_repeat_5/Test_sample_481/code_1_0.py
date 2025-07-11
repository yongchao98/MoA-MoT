def find_different_imprints():
    """
    This function identifies and prints the different numbers on Teva brand
    fluoxetine 10mg and 20mg capsules.
    """
    # Imprints found on the capsules
    imprint_10mg = "TEVA 3107"
    imprint_20mg = "TEVA 3108"

    # Extract the numerical digits from each imprint string into sets
    numbers_10mg = {char for char in imprint_10mg if char.isdigit()}
    numbers_20mg = {char for char in imprint_20mg if char.isdigit()}

    # Find the numbers that are unique to the 10mg capsule
    unique_to_10mg = numbers_10mg.difference(numbers_20mg)

    # Find the numbers that are unique to the 20mg capsule
    unique_to_20mg = numbers_20mg.difference(numbers_10mg)

    # The problem asks to output each number in the final equation.
    # We will print the findings clearly.
    if unique_to_10mg:
        num_10mg = list(unique_to_10mg)[0]
        print(f"The Teva 10mg fluoxetine capsule is imprinted with the number {num_10mg} that is not on the 20mg version.")

    if unique_to_20mg:
        num_20mg = list(unique_to_20mg)[0]
        print(f"The Teva 20mg fluoxetine capsule is imprinted with the number {num_20mg} that is not on the 10mg version.")

find_different_imprints()
<<<7 and 8>>>