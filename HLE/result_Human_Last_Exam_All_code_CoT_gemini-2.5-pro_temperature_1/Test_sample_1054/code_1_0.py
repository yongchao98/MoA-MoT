def find_max_diameter_cutoff():
    """
    Finds the maximum Aortomesenteric diameter cutoff that has a sensitivity
    greater than 60% and a specificity greater than 80%.
    """
    # Data extracted from the provided ROC curve analysis table.
    # Each tuple represents: (Diameter_cutoff, Sensitivity, Specificity)
    data_points = [
        (15.1, 6.7, 100.0),
        (13.9, 13.3, 100.0),
        (12.8, 20.0, 100.0),
        (12.0, 26.7, 100.0),
        (11.4, 33.3, 100.0),
        (10.9, 40.0, 100.0),
        (10.4, 46.7, 100.0),
        (9.8, 53.3, 95.8),
        (9.5, 60.0, 91.7),
        (9.2, 66.7, 87.5),
        (8.5, 73.3, 83.3),
        (7.9, 80.0, 75.0),
        (7.5, 86.7, 66.7),
        (6.9, 93.3, 54.2),
        (6.2, 100.0, 41.7),
        (5.5, 100.0, 29.2),
        (4.7, 100.0, 16.7),
        (4.2, 100.0, 4.2)
    ]

    # A list to store diameter cutoffs that meet the criteria.
    valid_diameters = []

    # Iterate through each data point to find ones that match the criteria.
    for diameter, sensitivity, specificity in data_points:
        # Check if sensitivity is > 60% AND specificity is > 80%.
        if sensitivity > 60 and specificity > 80:
            valid_diameters.append(diameter)

    # Find the maximum diameter from the list of valid diameters.
    if valid_diameters:
        max_valid_diameter = max(valid_diameters)
        # As per the instructions, print the final number of the result.
        print(max_valid_diameter)
    else:
        print("No diameter found that meets the specified criteria.")

# Execute the function to find and print the result.
find_max_diameter_cutoff()