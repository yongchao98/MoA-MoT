def find_largest_sunk_ship():
    """
    This function provides information about the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """

    ship_name = "Champlain"
    # Gross Register Tonnage (GRT) is a measure of a ship's internal volume.
    # Standard displacement for warships is measured differently (in tons of water displaced).
    # For a passenger liner like the Champlain, GRT is the standard measure of size.
    displacement = 28124
    unit = "GRT (Gross Register Tonnage)"
    cause = "a mine laid by the German submarine U-65"
    date_sunk = "June 17, 1940"

    print("Identifying the largest French ship by displacement sunk by a U-boat before the 1940 armistice:")
    print("-" * 80)
    print(f"The ship was the passenger liner '{ship_name}'.")
    print(f"Its displacement was {displacement} {unit}.")
    print(f"It was sunk on {date_sunk} after striking {cause}.")
    print("-" * 80)
    print("Final answer:")
    print(f"Ship Name: {ship_name}")
    print(f"Displacement: {displacement} {unit}")


if __name__ == "__main__":
    find_largest_sunk_ship()