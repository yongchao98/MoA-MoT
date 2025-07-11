def azobenzene_table_function():
    """
    This script describes the functionally accurate behavior of a picnic table
    shaped like an azobenzene molecule during a day/night cycle.
    """
    trans_form_description = "a straight, elongated shape (the 'trans-isomer')"
    cis_form_description = "a bent, V-shaped configuration (the 'cis-isomer')"

    print("--- Functional Cycle of the Azobenzene Picnic Table ---")
    
    # --- Sunrise ---
    print("\nEvent 1: Sunrise")
    print("The sun rises, exposing the table to UV radiation.")
    print(f"ACTION: To be functionally accurate, the table must change its structure.")
    print(f"The table converts from {trans_form_description} to {cis_form_description}.")

    # --- Sunset ---
    print("\nEvent 2: Sunset")
    print("The sun sets, and the source of UV radiation is removed.")
    print(f"ACTION: The table must now revert to its most stable state.")
    print(f"Overnight, the table slowly relaxes from {cis_form_description} back to {trans_form_description}.")
    
    print("\n--- Cycle Complete ---")

# Run the simulation
azobenzene_table_function()
