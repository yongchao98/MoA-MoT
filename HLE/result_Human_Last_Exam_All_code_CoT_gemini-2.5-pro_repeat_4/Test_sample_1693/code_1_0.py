def find_largest_sunk_ship():
    """
    This function identifies the largest French ship sunk by a U-boat's action
    before the armistice of June 22, 1940.
    """
    # Data format: (Ship Name, Type, Gross Register Tonnage)
    # This data represents a selection of the largest French ships sunk by U-boats
    # or mines laid by U-boats in the specified period.
    ships_data = [
        ("Émile Miguet", "Tanker", 14115),
        ("Cheyenne", "Tanker", 8825),
        ("Brazza", "Passenger Ship", 10387),
        ("Champlain", "Passenger Liner", 28124)
    ]

    # Initialize variables to store the details of the largest ship found.
    largest_ship_name = ""
    largest_ship_displacement = 0

    # Iterate through the data to find the ship with the maximum displacement.
    for name, ship_type, displacement in ships_data:
        if displacement > largest_ship_displacement:
            largest_ship_displacement = displacement
            largest_ship_name = name

    # Output the result.
    if largest_ship_name:
        # The final answer includes the name and the key number (displacement).
        print(f"The largest French ship sunk by a U-boat's action before the 1940 armistice was the S.S. {largest_ship_name}.")
        # The number from the final 'equation' or finding is the displacement.
        print(f"Its displacement was {largest_ship_displacement} tons.")
    else:
        print("Could not determine the largest ship from the provided data.")

# Execute the function to find and print the answer.
find_largest_sunk_ship()