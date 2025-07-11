def find_imprint_differences():
    """
    Identifies and prints the differing numbers on two Teva fluoxetine capsules.
    """
    # Numerical imprints for Teva fluoxetine 10mg and 20mg capsules
    imprint_10mg = "3105"
    imprint_20mg = "3106"

    diff_10mg_list = []
    diff_20mg_list = []

    # Compare the imprints digit by digit
    # Assuming imprints are the same length
    for i in range(len(imprint_10mg)):
        if imprint_10mg[i] != imprint_20mg[i]:
            diff_10mg_list.append(imprint_10mg[i])
            diff_20mg_list.append(imprint_20mg[i])

    # Format the imprints with spaces for the final output "equation"
    spaced_imprint_10mg = " ".join(imprint_10mg)
    spaced_imprint_20mg = " ".join(imprint_20mg)

    print(f"The number on the Teva 10mg fluoxetine capsule is: {spaced_imprint_10mg}")
    print(f"The number on the Teva 20mg fluoxetine capsule is: {spaced_imprint_20mg}")
    
    if diff_10mg_list:
        print("\nThe different numbers are:")
        print(f"From the 10mg capsule: {', '.join(diff_10mg_list)}")
        print(f"From the 20mg capsule: {', '.join(diff_20mg_list)}")
    else:
        print("\nThe numerical imprints are identical.")

find_imprint_differences()
