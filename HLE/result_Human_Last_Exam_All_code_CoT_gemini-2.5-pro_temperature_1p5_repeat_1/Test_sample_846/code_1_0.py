def simulate_azobenzene_table():
    """
    Simulates the daily functional behavior of an azobenzene picnic table.
    """

    # Define the properties of the two isomer states
    trans_isomer_shape = "flat (trans-isomer)"
    cis_isomer_shape = "bent (cis-isomer)"

    # Define the wavelengths that trigger the transformations
    uv_light_wavelength = 365  # Wavelength in nanometers for trans -> cis
    visible_light_wavelength = 440 # Wavelength in nanometers for cis -> trans

    # --- Sunrise Simulation ---
    print("--- AT SUNRISE ---")
    print("The sun appears, exposing the table to UV light.")
    print("The table must change from its stable, flat shape to a bent shape.")
    print("This transformation is called photoisomerization.")
    print("\nThe functional equation for this process is:")
    # Print the equation with the specific number for the wavelength
    print(f"UV Light ({uv_light_wavelength} nm) + {trans_isomer_shape} -> {cis_isomer_shape}")

    print("\n" + "="*50 + "\n")

    # --- Sunset Simulation ---
    print("--- AT SUNSET ---")
    print("The sun disappears, removing the UV light source.")
    print("The table must return to its stable, flat shape.")
    print("This transformation is called thermal relaxation (and can be sped up by visible light).")
    print("\nThe functional equation for this process is:")
    # Print the equation representing the return to the stable state
    print(f"Darkness/Heat + {cis_isomer_shape} -> {trans_isomer_shape}")
    print(f"(Alternatively: Visible Light ({visible_light_wavelength} nm) + {cis_isomer_shape} -> {trans_isomer_shape})")


# Run the simulation
simulate_azobenzene_table()