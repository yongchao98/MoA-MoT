def find_imprint_differences():
    """
    This script identifies and compares the numerical imprints on two
    different strengths of Teva brand fluoxetine capsules.
    """
    # Based on drug identification data, these are the numerical imprints.
    # The imprint on Teva 10mg fluoxetine capsules is "TEVA 3108".
    imprint_10mg = "3108"

    # The imprint on Teva 20mg fluoxetine capsules is "TEVA 3109".
    imprint_20mg = "3109"

    # Lists to store the digits that are different
    diff_digits_from_10mg = []
    diff_digits_from_20mg = []

    # Iterate through the strings and compare them character by character
    # This assumes the imprint numbers are the same length.
    for i in range(len(imprint_10mg)):
        if imprint_10mg[i] != imprint_20mg[i]:
            diff_digits_from_10mg.append(imprint_10mg[i])
            diff_digits_from_20mg.append(imprint_20mg[i])
    
    # Format the final output
    print(f"The number imprinted on the 10mg capsule is {imprint_10mg}.")
    print(f"The number imprinted on the 20mg capsule is {imprint_20mg}.")

    if diff_digits_from_10mg:
        # Construct the final "equation" or comparison statement
        diff1 = diff_digits_from_10mg[0]
        diff2 = diff_digits_from_20mg[0]
        print(f"Comparing {imprint_10mg} and {imprint_20mg}, the different numbers are {diff1} and {diff2}.")
    else:
        print("The numbers are identical.")

find_imprint_differences()