def simulate_azobenzene_table():
    """
    This script describes the functionally accurate behavior of a picnic table
    shaped like an azobenzene molecule throughout a day.
    """
    # Define the two molecular states and their corresponding table shapes
    trans_state = "stable, flat, extended shape (trans-isomer)"
    cis_state = "less-stable, bent, V-shape (cis-isomer)"

    print("--- Simulating a day for the Azobenzene Picnic Table ---")

    # Sunrise event
    print("\nEvent: SUNRISE")
    print("Action: The sun's UV light hits the table.")
    print(f"The table changes from its {trans_state}...")
    print(f"...to its {cis_state}.")
    print("Result: The table bends or contracts.")

    # Sunset event
    print("\nEvent: SUNSET")
    print("Action: UV light fades and darkness falls.")
    print(f"The table thermally relaxes from its {cis_state}...")
    print(f"...back to its {trans_state}.")
    print("Result: The table straightens out or extends.")

if __name__ == '__main__':
    simulate_azobenzene_table()