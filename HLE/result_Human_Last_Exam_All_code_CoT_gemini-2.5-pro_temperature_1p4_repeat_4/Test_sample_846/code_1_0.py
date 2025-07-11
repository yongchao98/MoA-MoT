def simulate_azobenzene_table():
    """
    This function simulates the shape-shifting behavior of a picnic table
    designed to mimic the photoisomerization of an azobenzene molecule.
    """
    # The 'trans' isomer is more stable and is the default state in the dark.
    # We will represent this as the table's "flat" state.
    table_state = "trans (flat)"
    print(f"Pre-dawn: In the dark, the table is in its most stable state: {table_state}")
    print("-" * 50)

    # At sunrise, sunlight provides UV radiation.
    print("Event: The sun rises, exposing the table to UV light.")
    # UV light triggers the isomerization from the 'trans' to the 'cis' state.
    # We will represent the bent 'cis' state as the table's "bent" state.
    table_state = "cis (bent)"
    print(f"Result: The table bends and changes its conformation.")
    print(f"Current state: {table_state}")
    print("-" * 50)

    # At sunset, the UV source is removed, and the molecule thermally relaxes.
    print("Event: The sun sets, and darkness falls.")
    # The 'cis' isomer reverts to the more stable 'trans' isomer.
    table_state = "trans (flat)"
    print(f"Result: The table straightens out, returning to its stable form.")
    print(f"Final state: {table_state}")

simulate_azobenzene_table()